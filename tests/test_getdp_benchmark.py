import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
GETDP = os.environ.get("GETDP_EXECUTABLE") or shutil.which("getdp")
if not GETDP and (ROOT / ".tools/getdp-git-MacOSARM/bin/getdp").is_file():
    GETDP = str(ROOT / ".tools/getdp-git-MacOSARM/bin/getdp")
AVAILABLE = GETDP and all(importlib.util.find_spec(p) for p in ("dolfinx", "gmsh", "scipy"))
if AVAILABLE:
    import numpy as np
    from dolfinx import mesh as dmesh
    from dolfinx.io import XDMFFile, gmsh as gmshio
    from mpi4py import MPI
    from examples.fem_tank import solve_cem
    from examples.getdp_cem import inspect_mesh, solve_getdp
    from examples.getdp_benchmark import write_mesh, compare_solvers, write_benchmark


@unittest.skipUnless(AVAILABLE, "GetDP plus DOLFINx/Gmsh environment required")
class IndependentCEMTests(unittest.TestCase):
    def test_same_mesh_scalar_and_off_diagonal_tensor_p1_p2(self):
        sigma_tensor = np.array([[.5, .12, -.05], [.12, .8, .09], [-.05, .09, .4]])
        for degree in (1, 2):
            for sigma in (.5, sigma_tensor):
                with self.subTest(degree=degree, sigma=sigma):
                    row = compare_solvers(ROOT / "data/synthetic_fem/finite_patch_0.msh", sigma, degree=degree)
                    for metric in ("voltage_relative_difference", "potential_sample_relative_rms",
                                   "electric_sample_relative_rms", "current_sample_relative_rms",
                                   "bulk_power_relative_difference", "contact_power_relative_difference",
                                   "getdp_relative_power_error"):
                        self.assertLess(row[metric], 1e-8, metric)
                    self.assertLess(row["getdp_max_current_error_a"], 1e-11)

    def test_adaptive_geometry_connectivity_tags_and_current_sign_round_trip(self):
        with XDMFFile(MPI.COMM_SELF, str(ROOT / "data/synthetic_fem_adaptive/mesh.xdmf"), "r") as file:
            domain = file.read_mesh(name="mesh")
            domain.topology.create_connectivity(2, 3)
            data = SimpleNamespace(mesh=domain, cell_tags=file.read_meshtags(domain, name="materials"),
                                   facet_tags=file.read_meshtags(domain, name="boundaries"))

        def geometry_signature(data, tags):
            vertices = dmesh.entities_to_geometry(data.mesh, tags.dim, tags.indices)
            # Ignore local numbering and orientation; preserve every geometric vertex and tag.
            return sorted((int(tag), tuple(sorted(tuple(x) for x in data.mesh.geometry.x[v])))
                          for v, tag in zip(vertices, tags.values))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "adaptive.msh"
            write_mesh(data, path)
            loaded = gmshio.read_from_msh(str(path), MPI.COMM_SELF)
            self.assertEqual(geometry_signature(data, data.cell_tags), geometry_signature(loaded, loaded.cell_tags))
            self.assertEqual(geometry_signature(data, data.facet_tags), geometry_signature(loaded, loaded.facet_tags))
            materials, areas = inspect_mesh(path)
            self.assertEqual(materials, {1})
            for tag in (101, 102):
                self.assertAlmostEqual(areas[tag], .0012, places=14)
            currents, contacts = {101: .001, 102: -.001}, {101: .002, 102: .005}
            original = solve_cem(data, .5, currents, contacts, degree=2)
            independent = solve_getdp(path, .5, currents, contacts, degree=2, executable=GETDP)
            self.assertAlmostEqual(independent.electrode_potential_v[101], original.electrode_potential_v[101], places=10)
            for tag, current in currents.items():
                self.assertAlmostEqual(independent.electrode_current_a[tag], current, places=11)

    def test_release_hashes_native_evidence_and_energy_interval(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = write_benchmark(root)
            self.assertEqual(len(rows), 9)
            final = next(row for row in rows if row["case"] == "adaptive")
            resistance = final["getdp_voltage_v"] / .001
            self.assertGreaterEqual(resistance, final["resistance_lower_ohm"] - 1e-8)
            self.assertLessEqual(resistance, final["resistance_upper_ohm"] + 1e-8)
            refined = [row for row in rows if row["case"].startswith("finite_patch") and row["degree"] == 2]
            self.assertEqual(len(refined), 3)
            self.assertLess(abs(refined[-1]["getdp_voltage_v"] - final["getdp_voltage_v"]),
                            abs(refined[0]["getdp_voltage_v"] - final["getdp_voltage_v"]))
            for row in rows:
                self.assertLess(row["current_sample_relative_rms"], 1e-8)
            manifest = json.loads((root / "manifest.json").read_text())
            for path, digest in manifest["files_sha256"].items():
                self.assertEqual(hashlib.sha256((root / path).read_bytes()).hexdigest(), digest)
            for name in ("problem.res", "problem.pro", "parameters.pro", "solver-info.txt", "solver.log"):
                self.assertTrue((root / "adaptive_p2" / name).is_file())
            with np.load(root / "adaptive_p2/comparison.npz", allow_pickle=False) as fields:
                self.assertEqual(fields["getdp_j_a_m2"].shape, (6826, 3))
                np.testing.assert_allclose(fields["getdp_j_a_m2"], fields["dolfinx_j_a_m2"], atol=1e-8)
