"""Physical oracles for an independent tetrahedral complete-electrode model."""

import importlib.util
import tempfile
from pathlib import Path
import unittest

AVAILABLE = all(importlib.util.find_spec(p) is not None
                for p in ("dolfinx", "gmsh", "scipy", "mpi4py"))
if AVAILABLE:
    import numpy as np
    from dolfinx import fem, mesh
    from mpi4py import MPI
    import ufl
    from examples.fem_tank import build_tank_mesh, solve_cem


@unittest.skipUnless(AVAILABLE, "Optional DOLFINx/Gmsh environment required")
class CompleteElectrodeTests(unittest.TestCase):
    def setUp(self):
        self.lengths = np.array([.12, .08, .06])
        self.area = .08 * .06
        self.currents = {101: .001, 102: -.001}  # inward-positive
        self.contacts = {101: .002, 102: .005}  # ohm m^2, not ohm

    def test_cad_tags_volume_area_and_round_trip(self):
        from dolfinx.io import gmsh as gmshio
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tank.msh"
            data = build_tank_mesh(patch_fraction=.5, layered=True, mesh_path=path)
            self.assertTrue(path.is_file())
            restored = gmshio.read_from_msh(str(path), MPI.COMM_SELF)
            for item in (data, restored):
                dx = ufl.Measure("dx", domain=item.mesh, subdomain_data=item.cell_tags)
                ds = ufl.Measure("ds", domain=item.mesh, subdomain_data=item.facet_tags)
                for tag in (1, 2):
                    self.assertAlmostEqual(fem.assemble_scalar(fem.form(1 * dx(tag))),
                                           np.prod(self.lengths) / 2, places=13)
                for tag in (101, 102):
                    self.assertAlmostEqual(fem.assemble_scalar(fem.form(1 * ds(tag))),
                                           self.area * .5**2, places=13)
                self.assertEqual(set(item.facet_tags.values), {101, 102, 199})

    def test_slab_analytic_voltage_contact_power_and_fem_order(self):
        data = build_tank_mesh()
        j = .001 / self.area
        resistance = (.12 / .5 + .002 + .005) / self.area
        for degree in (1, 2):
            with self.subTest(degree=degree):
                out = solve_cem(data, .5, self.currents, self.contacts, degree=degree)
                self.assertAlmostEqual(out.electrode_potential_v[101], .001 * resistance, places=11)
                self.assertAlmostEqual(out.electrode_potential_v[102], 0., places=14)
                x = out.potential_v.function_space.tabulate_dof_coordinates()[:, 0]
                expected = j * ((.12 - x) / .5 + .005)
                np.testing.assert_allclose(out.potential_v.x.array, expected, atol=2e-11)
                self.assertAlmostEqual(out.bulk_power_w, .001**2 * .12 / (.5 * self.area), places=13)
                self.assertAlmostEqual(out.contact_power_w, .001**2 * (.002 + .005) / self.area, places=13)
                self.check_balance(out)

    def test_layered_slab_matches_series_resistance_and_fv(self):
        from examples.conductive_tank import solve_plate_current
        data = build_tank_mesh(layered=True)
        out = solve_cem(data, {1: .2, 2: .8}, self.currents, self.contacts)
        sigma = np.ones((8, 6, 4)) * .8
        sigma[:4] = .2
        fv = solve_plate_current(sigma, self.lengths / [8, 6, 4], .001)
        contact_drop = .001 * sum(self.contacts.values()) / self.area
        expected = .001 * (.06 / .2 + .06 / .8) / self.area
        self.assertAlmostEqual(out.electrode_potential_v[101] - contact_drop, expected, places=11)
        self.assertAlmostEqual(out.electrode_potential_v[101] - contact_drop,
                               float(fv.boundary_potential_v["x-"][0, 0]), places=11)
        self.check_balance(out)

    def test_tensor_conductivity_and_rotation_covariance(self):
        data = build_tank_mesh()
        tensor = np.diag([.3, .6, .9])
        out = solve_cem(data, tensor, self.currents, self.contacts, degree=2)
        theta = .61
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                             [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
        data.mesh.geometry.x[:] = data.mesh.geometry.x @ rotation.T
        rotated = solve_cem(data, rotation @ tensor @ rotation.T,
                            self.currents, self.contacts, degree=2)
        np.testing.assert_allclose(rotated.potential_v.x.array, out.potential_v.x.array, atol=2e-11)
        expected_drop = .001 * (.12 / .3 + sum(self.contacts.values())) / self.area
        self.assertAlmostEqual(rotated.electrode_potential_v[101], expected_drop, places=11)
        j = rotated.current_density_a_m2.x.array.reshape(-1, 3)
        np.testing.assert_allclose(j, np.broadcast_to((.001 / self.area) * rotation[:, 0], j.shape),
                                   atol=2e-10)
        self.check_balance(rotated)

    def test_finite_patches_polarity_off_gauge_and_reciprocity(self):
        data = build_tank_mesh(patch_fraction=.5)
        # A third, floating equipotential electrode covers the remaining outer surface.
        contacts = self.contacts.copy()
        contacts[199] = .004
        drive_a = {101: .001, 102: -.001, 199: 0.}
        drive_b = {101: 0., 102: -.001, 199: .001}
        a = solve_cem(data, .5, drive_a, contacts)
        b = solve_cem(data, .5, drive_b, contacts)
        mutual_ab = sum(drive_a[t] * b.electrode_potential_v[t] for t in drive_a)
        mutual_ba = sum(drive_b[t] * a.electrode_potential_v[t] for t in drive_a)
        self.assertAlmostEqual(mutual_ab, mutual_ba, places=13)
        for out in (a, b):
            self.check_balance(out)
        reverse = solve_cem(data, .5, {t: -v for t, v in drive_a.items()}, contacts)
        np.testing.assert_allclose(reverse.potential_v.x.array, -a.potential_v.x.array, atol=1e-12)
        zero = solve_cem(data, .5, dict.fromkeys(drive_a, 0.), contacts)
        np.testing.assert_array_equal(zero.potential_v.x.array, 0.)
        shifted = solve_cem(data, .5, drive_a, contacts, reference_tag=101)
        np.testing.assert_allclose(shifted.potential_v.x.array,
                                   a.potential_v.x.array - a.electrode_potential_v[101], atol=2e-11)

    def test_invalid_conductivity_contact_current_and_tags(self):
        data = build_tank_mesh()
        bad = [(-.1, self.currents, self.contacts),
               (np.diag([1., -1., 1.]), self.currents, self.contacts),
               (np.array([[1., .3, 0], [0, 1., 0], [0, 0, 1.]]), self.currents, self.contacts),
               (complex(.5, .1), self.currents, self.contacts),
               (.5, {101: .001, 102: -.0009}, self.contacts),
               (.5, {101: .001 + .1j, 102: -.001 - .1j}, self.contacts),
               (.5, self.currents, {101: 0., 102: .002}),
               (.5, self.currents, {101: np.nan, 102: .002}),
               (.5, self.currents, {101: .002}),
               ({2: .5}, self.currents, self.contacts),
               (.5, {777: .001, 102: -.001}, {777: .002, 102: .002})]
        for sigma, currents, contacts in bad:
            with self.subTest(sigma=sigma, currents=currents, contacts=contacts):
                with self.assertRaises(ValueError):
                    solve_cem(data, sigma, currents, contacts)
        for degree in (0, 1.5):
            with self.assertRaises(ValueError):
                solve_cem(data, .5, self.currents, self.contacts, degree=degree)
        with self.assertRaises(ValueError):
            solve_cem(data, .5, self.currents, self.contacts, reference_tag=999)

    def test_invalid_geometry(self):
        for kwargs in ({"lengths_m": [0, 1, 1]}, {"mesh_size_m": -1},
                       {"patch_fraction": 0}, {"patch_fraction": 1.1},
                       {"lengths_m": [1, 2]}, {"mesh_size_m": np.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                build_tank_mesh(**kwargs)

    def check_balance(self, out):
        for tag, applied in out.applied_current_a.items():
            self.assertAlmostEqual(out.electrode_current_a[tag], applied, places=12)
        self.assertLess(out.relative_residual, 1e-10)
        self.assertAlmostEqual(sum(out.electrode_current_a.values()), 0., places=12)
        self.assertAlmostEqual(out.supplied_power_w, out.bulk_power_w + out.contact_power_w, places=12)


if __name__ == "__main__":
    unittest.main()
