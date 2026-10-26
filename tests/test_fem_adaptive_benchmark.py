import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

AVAILABLE = all(importlib.util.find_spec(p) is not None for p in ("dolfinx", "gmsh", "scipy"))
if AVAILABLE:
    import numpy as np
    from dolfinx.io import XDMFFile
    from mpi4py import MPI
    from examples.fem_tank import solve_cem
    from examples.fem_adaptive_benchmark import write_benchmark


@unittest.skipUnless(AVAILABLE, "Optional DOLFINx/Gmsh environment required")
class AdaptiveReleaseTests(unittest.TestCase):
    def test_release_fields_hashes_and_mesh_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_benchmark(root)
            parameters = json.loads((root / "parameters.json").read_text())
            self.assertTrue(parameters["converged"])
            self.assertEqual(parameters["stop_reason"], "tolerance")
            with (root / "history.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertGreaterEqual(len(rows), 2)
            self.assertLess(float(rows[-1]["relative_energy_gap"]), .01)
            with np.load(root / "fields.npz", allow_pickle=False) as data:
                cells = int(rows[-1]["cells"])
                self.assertEqual(data["rt_cell_dofs"].shape, (cells, 15))
                self.assertEqual(data["cell_permutations"].shape, (cells,))
                self.assertEqual(data["cell_center_j_a_m2"].shape, (cells, 3))
                self.assertAlmostEqual(sum(data["cell_indicator_squared_w"]), float(rows[-1]["estimator_squared_w"]), places=12)
                self.assertLess(np.max(abs(data["cell_current_imbalance_a"])), 1e-11)
            manifest = json.loads((root / "manifest.json").read_text())
            for path, digest in manifest["files_sha256"].items():
                self.assertEqual(hashlib.sha256((root / path).read_bytes()).hexdigest(), digest)
            with XDMFFile(MPI.COMM_SELF, str(root / "mesh.xdmf"), "r") as file:
                domain = file.read_mesh(name="mesh")
                domain.topology.create_connectivity(2, 3)
                materials = file.read_meshtags(domain, name="materials")
                boundaries = file.read_meshtags(domain, name="boundaries")
            restored = SimpleNamespace(mesh=domain, cell_tags=materials, facet_tags=boundaries)
            solution = solve_cem(restored, .5, {101: .001, 102: -.001}, {101: .002, 102: .005}, degree=2)
            self.assertAlmostEqual(solution.supplied_power_w, float(rows[-1]["primal_power_w"]), places=12)

    def test_limited_release_reports_unmet_tolerance_and_repeats(self):
        with tempfile.TemporaryDirectory() as directory:
            roots = [Path(directory) / name for name in ("a", "b")]
            for root in roots:
                write_benchmark(root, max_refinements=0)
                parameters = json.loads((root / "parameters.json").read_text())
                self.assertFalse(parameters["converged"])
                self.assertEqual(parameters["stop_reason"], "refinement_limit")
            self.assertEqual((roots[0] / "history.csv").read_bytes(), (roots[1] / "history.csv").read_bytes())
            with np.load(roots[0] / "fields.npz") as a, np.load(roots[1] / "fields.npz") as b:
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key])
