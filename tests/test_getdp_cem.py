import importlib.util
import os
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
GETDP = os.environ.get("GETDP_EXECUTABLE") or shutil.which("getdp")
LOCAL = ROOT / ".tools/getdp-git-MacOSARM/bin/getdp"
if not GETDP and LOCAL.is_file():
    GETDP = str(LOCAL)
AVAILABLE = GETDP and all(importlib.util.find_spec(p) is not None for p in ("numpy", "gmsh"))
if AVAILABLE:
    import numpy as np
    from examples.getdp_cem import solve_getdp


@unittest.skipUnless(AVAILABLE, "GetDP executable plus NumPy/Gmsh required")
class GetDPCEMTests(unittest.TestCase):
    def setUp(self):
        self.mesh = ROOT / "data/synthetic_fem/homogeneous.msh"
        self.currents = {101: .001, 102: -.001}
        self.contacts = {101: .002, 102: .005}

    def solve(self, sigma=.5, **kwargs):
        return solve_getdp(self.mesh, sigma, self.currents, self.contacts, executable=GETDP, **kwargs)

    def check_balance(self, out):
        for tag, applied in out.applied_current_a.items():
            self.assertAlmostEqual(out.electrode_current_a[tag], applied, places=11)
        self.assertAlmostEqual(out.bulk_power_w + out.contact_power_w, out.supplied_power_w, places=11)

    def test_analytic_slab_fields_voltage_and_contact_power_p1_p2(self):
        area = .08 * .06
        j = .001 / area
        for degree in (1, 2):
            with self.subTest(degree=degree):
                out = self.solve(degree=degree)
                self.assertAlmostEqual(out.electrode_potential_v[101], j * (.12 / .5 + .007), places=10)
                self.assertAlmostEqual(out.electrode_potential_v[102], 0., places=13)
                expected = j * ((.12 - out.sample_coordinates_m[:, 0]) / .5 + .005)
                np.testing.assert_allclose(out.potential_v, expected, atol=1e-10)
                np.testing.assert_allclose(out.current_density_a_m2,
                                           np.broadcast_to([j, 0., 0.], out.current_density_a_m2.shape), atol=2e-9)
                self.assertAlmostEqual(out.bulk_power_w, .001**2 * .12 / (.5 * area), places=12)
                self.assertAlmostEqual(out.contact_power_w, .001**2 * .007 / area, places=12)
                self.check_balance(out)

    def test_layered_series_resistance(self):
        self.mesh = ROOT / "data/synthetic_fem/layered.msh"
        out = self.solve({1: .2, 2: .8}, degree=2)
        self.assertAlmostEqual(out.electrode_potential_v[101], .001 * (.06 / .2 + .06 / .8 + .007) / (.08 * .06), places=10)
        self.check_balance(out)

    def test_diagonal_tensor_slab(self):
        out = self.solve(np.diag([.3, .6, .9]), degree=2)
        self.assertAlmostEqual(out.electrode_potential_v[101], .001 * (.12 / .3 + .007) / (.08 * .06), places=10)
        self.check_balance(out)

    def test_floating_electrode_reciprocity_polarity_gauge_and_off(self):
        self.mesh = ROOT / "data/synthetic_fem/finite_patch_0.msh"
        self.currents[199] = 0.
        self.contacts[199] = .004
        a = self.solve()
        b_currents = {101: 0., 102: -.001, 199: .001}
        b = solve_getdp(self.mesh, .5, b_currents, self.contacts, executable=GETDP)
        ab = sum(self.currents[t] * b.electrode_potential_v[t] for t in self.currents)
        ba = sum(b_currents[t] * a.electrode_potential_v[t] for t in self.currents)
        self.assertAlmostEqual(ab, ba, places=11)
        self.check_balance(a)
        self.check_balance(b)
        for scale, reference in ((-1., 199), (0., 199), (1., 101)):
            out = solve_getdp(self.mesh, .5, {t: scale * i for t, i in self.currents.items()},
                              self.contacts, reference_tag=reference, executable=GETDP)
            shift = scale * a.electrode_potential_v[reference]
            np.testing.assert_allclose(out.potential_v, scale * a.potential_v - shift, atol=1e-9)
            np.testing.assert_allclose(out.current_density_a_m2, scale * a.current_density_a_m2, atol=1e-8)
            self.check_balance(out)

    def test_native_problem_outputs_and_provenance_are_retained(self):
        with tempfile.TemporaryDirectory() as directory:
            out = self.solve(output_dir=directory)
            for name in ("problem.pro", "problem.res", "solver.log", "solver-info.txt"):
                self.assertTrue((Path(directory) / name).is_file())
            self.assertIn("GetDP", out.solver_info)
            self.assertGreater(len(out.sample_coordinates_m), 10)

    def test_invalid_physics_order_and_mesh_are_rejected(self):
        for sigma in (-.5, np.diag([1., -1., 1.]), complex(.5, .1),
                      np.array([[1, .2, 0], [0, 1, 0], [0, 0, 1]]), {9: .5}):
            with self.subTest(sigma=sigma), self.assertRaises(ValueError):
                self.solve(sigma)
        for degree in (0, 3, 1.5, True):
            with self.subTest(degree=degree), self.assertRaises(ValueError):
                self.solve(degree=degree)
        for currents, contacts in (({101: .001, 102: -.0009}, self.contacts),
                                    (self.currents, {101: 0., 102: .005}),
                                    (self.currents, {101: .002}),
                                    ({777: .001, 102: -.001}, {777: .002, 102: .005})):
            with self.assertRaises(ValueError):
                solve_getdp(self.mesh, .5, currents, contacts, executable=GETDP)
        with self.assertRaises(ValueError):
            self.solve(reference_tag=777)
