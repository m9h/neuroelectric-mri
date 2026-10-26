"""Analytical and independent weak/strong checks of equilibrated FEM currents."""

import importlib.util
import unittest

AVAILABLE = all(importlib.util.find_spec(p) is not None for p in ("dolfinx", "gmsh", "scipy"))
if AVAILABLE:
    import numpy as np
    import ufl
    from dolfinx import fem
    from examples.fem_tank import build_tank_mesh, solve_cem
    from examples.fem_flux import reconstruct_current, mark_cells


@unittest.skipUnless(AVAILABLE, "Optional DOLFINx/Gmsh environment required")
class ConservativeCurrentTests(unittest.TestCase):
    def setUp(self):
        self.currents = {101: .001, 102: -.001}
        self.contacts = {101: .002, 102: .005}

    def solve(self, data, sigma=.5, **kwargs):
        return solve_cem(data, sigma, self.currents, self.contacts, **kwargs)

    def check_conservation(self, data, flux):
        j = flux.current_density_a_m2
        n = ufl.FacetNormal(data.mesh)
        dx = ufl.Measure("dx", domain=data.mesh)
        dS = ufl.Measure("dS", domain=data.mesh)
        ds = ufl.Measure("ds", domain=data.mesh, subdomain_data=data.facet_tags)
        # Check the returned finite element function independently of its diagnostics.
        divergence_l2 = np.sqrt(abs(fem.assemble_scalar(fem.form(ufl.div(j)**2 * dx))))
        jump_l2 = np.sqrt(abs(fem.assemble_scalar(fem.form(ufl.jump(j, n)**2 * dS))))
        self.assertLess(divergence_l2, 2e-9)
        self.assertLess(jump_l2, 2e-12)
        for tag, applied in self.currents.items():
            inward = -fem.assemble_scalar(fem.form(ufl.dot(j, n) * ds(tag)))
            self.assertAlmostEqual(inward, applied, places=11)
        self.assertLess(np.max(abs(flux.cell_current_imbalance_a)), 1e-11)
        self.assertLess(flux.relative_constraint_residual, 1e-9)

    def test_slab_exact_current_and_zero_estimator_both_rt_orders(self):
        data = build_tank_mesh()
        solution = self.solve(data, degree=2)
        for degree in (1, 2):
            with self.subTest(degree=degree):
                flux = reconstruct_current(data, solution, degree=degree)
                j = flux.current_density_a_m2
                exact = ufl.as_vector([.001 / (.08 * .06), 0., 0.])
                error = fem.assemble_scalar(fem.form(ufl.inner(j - exact, j - exact) * ufl.dx))
                self.assertLess(abs(error), 1e-22)
                self.assertLess(flux.estimator_squared_w, 1e-20)
                self.assertAlmostEqual(flux.complementary_power_w, solution.supplied_power_w, places=12)
                self.check_conservation(data, flux)

    def test_layered_rotated_tensor_has_continuous_correct_current(self):
        data = build_tank_mesh(layered=True)
        theta = .43
        r = np.array([[np.cos(theta), -np.sin(theta), 0],
                      [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
        data.mesh.geometry.x[:] = data.mesh.geometry.x @ r.T
        sigma = {1: r @ np.diag([.2, .4, .6]) @ r.T,
                 2: r @ np.diag([.8, .5, .3]) @ r.T}
        solution = self.solve(data, sigma, degree=2)
        flux = reconstruct_current(data, solution, degree=2)
        exact = ufl.as_vector((.001 / (.08 * .06) * r[:, 0]).tolist())
        error = fem.assemble_scalar(fem.form(ufl.inner(flux.current_density_a_m2 - exact,
                                                       flux.current_density_a_m2 - exact) * ufl.dx))
        self.assertLess(abs(error), 1e-22)
        self.check_conservation(data, flux)

    def test_finite_patches_remove_jumps_and_satisfy_contact_energy_identity(self):
        data = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5, layered=True)
        solution = self.solve(data, {1: .2, 2: .8}, degree=2)
        n = ufl.FacetNormal(data.mesh)
        raw_jump = fem.assemble_scalar(fem.form(ufl.jump(solution.current_density_a_m2, n)**2 * ufl.dS))
        self.assertGreater(raw_jump, 1e-8)  # the pre-existing gradient field is not equilibrated
        for degree in (1, 2):
            with self.subTest(degree=degree):
                flux = reconstruct_current(data, solution, degree=degree)
                self.check_conservation(data, flux)
                ds = ufl.Measure("ds", domain=data.mesh, subdomain_data=data.facet_tags)
                leakage = fem.assemble_scalar(fem.form(ufl.dot(flux.current_density_a_m2, n)**2 * ds(199)))
                self.assertLess(abs(leakage), 1e-22)
                self.assertGreater(flux.bulk_mismatch_w, 0.)
                self.assertGreater(flux.contact_mismatch_w, 0.)
                self.assertAlmostEqual(flux.estimator_squared_w,
                                       flux.bulk_mismatch_w + flux.contact_mismatch_w, places=13)
                self.assertAlmostEqual(flux.complementary_power_w - solution.supplied_power_w,
                                       flux.estimator_squared_w, places=12)
                self.assertAlmostEqual(sum(flux.cell_indicator_squared_w), flux.estimator_squared_w, places=12)
                self.assertTrue(np.all(flux.cell_indicator_squared_w >= 0))

    def test_energy_interval_contains_independent_finer_primal_result(self):
        coarse = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5)
        solution = self.solve(coarse, degree=1)
        flux = reconstruct_current(coarse, solution, degree=2)
        fine = build_tank_mesh(mesh_size_m=.02, patch_fraction=.5)
        reference = self.solve(fine, degree=2)
        self.assertLess(solution.supplied_power_w, reference.supplied_power_w)
        self.assertGreater(flux.complementary_power_w, reference.supplied_power_w)
        self.assertLess(reference.supplied_power_w - solution.supplied_power_w, flux.estimator_squared_w)

    def test_polarity_gauge_off_and_floating_electrodes(self):
        data = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5)
        contacts = {**self.contacts, 199: .004}
        currents = {**self.currents, 199: 0.}
        solution = solve_cem(data, .5, currents, contacts)
        baseline = reconstruct_current(data, solution)
        for sign, reference_tag in ((-1., 199), (1., 101), (0., 199)):
            driven = solve_cem(data, .5, {t: sign * i for t, i in currents.items()},
                               contacts, reference_tag=reference_tag)
            flux = reconstruct_current(data, driven)
            np.testing.assert_allclose(flux.current_density_a_m2.x.array,
                                       sign * baseline.current_density_a_m2.x.array, atol=1e-11)
            self.assertAlmostEqual(flux.estimator_squared_w, sign**2 * baseline.estimator_squared_w, places=12)
            for tag, applied in driven.applied_current_a.items():
                self.assertAlmostEqual(flux.electrode_current_a[tag], applied, places=11)

    def test_reconstruction_rejects_wrong_mesh_and_order(self):
        data = build_tank_mesh()
        solution = self.solve(data)
        for degree in (0, 3, 1.5, True):
            with self.subTest(degree=degree), self.assertRaises(ValueError):
                reconstruct_current(data, solution, degree=degree)
        with self.assertRaises(ValueError):
            reconstruct_current(build_tank_mesh(), solution)

    def test_bulk_marking_is_minimal_deterministic_and_has_zero_null(self):
        indicators = np.array([1., 4., 2., 4., 0.])
        np.testing.assert_array_equal(mark_cells(indicators, .6), [1, 3])
        np.testing.assert_array_equal(mark_cells(indicators, 1.), [1, 3, 2, 0])
        np.testing.assert_array_equal(mark_cells(np.zeros(5)), [])
        for bad in ([1., -1.], [np.nan], [[1., 2.]], [1 + 1j]):
            with self.assertRaises(ValueError):
                mark_cells(bad)
        for fraction in (0, -1, 1.1, np.nan):
            with self.assertRaises(ValueError):
                mark_cells(indicators, fraction)


if __name__ == "__main__":
    unittest.main()
