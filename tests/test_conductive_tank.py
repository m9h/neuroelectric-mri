"""Independent physical references for a conservative scalar tank solver."""

import importlib.util
import unittest

AVAILABLE = all(importlib.util.find_spec(n) is not None for n in ("numpy", "scipy"))
if AVAILABLE:
    import numpy as np
    from examples import conductive_tank as tank


@unittest.skipUnless(AVAILABLE, "Install requirements-tank.txt for the numerical tank")
class ConductiveTankTests(unittest.TestCase):
    def test_uniform_plate_solution_on_each_axis(self):
        shape = (6, 4, 3)
        spacing = np.array([.02, .015, .01])
        sigma, voltage = .5, .02
        for axis, name in enumerate("xyz"):
            with self.subTest(axis=name):
                result = tank.solve_tank(np.full(shape, sigma), spacing,
                                         {name + "-": voltage, name + "+": 0.})
                x = (np.arange(shape[axis]) + .5) * spacing[axis]
                length = shape[axis] * spacing[axis]
                expected = voltage * (1 - x / length)
                reshape = [1, 1, 1]
                reshape[axis] = shape[axis]
                np.testing.assert_allclose(result.potential_v,
                    np.broadcast_to(expected.reshape(reshape), shape), atol=1e-14)
                field = np.zeros((*shape, 3))
                field[..., axis] = voltage / length
                np.testing.assert_allclose(result.electric_field_v_m, field, atol=1e-12)
                np.testing.assert_allclose(result.current_density_a_m2, sigma * field, atol=1e-12)
                area = np.prod(np.delete(np.array(shape) * spacing, axis))
                current = sigma * voltage / length * area
                self.assertAlmostEqual(result.boundary_current_a[name + "+"], current, places=13)
                self.assertAlmostEqual(result.boundary_current_a[name + "-"], -current, places=13)
                for other in "xyz".replace(name, ""):
                    self.assertEqual(result.boundary_current_a[other + "-"], 0.)
                    self.assertEqual(result.boundary_current_a[other + "+"], 0.)
                self.assertAlmostEqual(result.dissipated_power_w, current * voltage, places=14)
                self.assertAlmostEqual(result.boundary_power_w, current * voltage, places=14)

    def test_one_cell_tank_has_two_half_cell_boundary_resistances(self):
        result = tank.solve_tank(np.full((1, 1, 1), .4), (.1, .02, .03),
                                 {"x-": .03, "x+": -.01})
        self.assertAlmostEqual(result.potential_v.item(), .01)
        self.assertAlmostEqual(result.boundary_current_a["x+"], .4 * .04 / .1 * .02 * .03)

    def test_layered_series_resistance_and_interface_flux(self):
        shape, dx = (10, 3, 2), .01
        sigma = np.full(shape, .2)
        sigma[5:] = .8
        voltage, area = .015, .03 * .02
        current_density = voltage / (.05 / .2 + .05 / .8)
        x = (np.arange(shape[0]) + .5) * dx
        resistance_to_x = np.minimum(x, .05) / .2 + np.maximum(x - .05, 0) / .8
        expected = voltage - current_density * resistance_to_x
        result = tank.solve_tank(sigma, (dx, dx, dx), {"x-": voltage, "x+": 0.})
        np.testing.assert_allclose(result.potential_v,
            np.broadcast_to(expected[:, None, None], shape), atol=1e-14)
        np.testing.assert_allclose(result.face_current_a[0], current_density * dx * dx, atol=1e-14)
        np.testing.assert_allclose(result.current_density_a_m2[..., 0], current_density, atol=1e-12)
        np.testing.assert_allclose(result.electric_field_v_m[..., 0], current_density / sigma, atol=1e-12)
        self.assertAlmostEqual(result.boundary_current_a["x+"], current_density * area, places=13)

    def test_cell_conservation_and_power_with_inclusion(self):
        sigma = np.full((7, 6, 5), .4)
        sigma[2:5, 2:4, 1:4] = .07
        spacing = (.01, .012, .009)
        original = sigma.copy()
        result = tank.solve_tank(sigma, spacing, {"x-": .02, "x+": 0.})
        # Compute conservation directly from exported face currents, not a solver residual.
        balance = sum(np.diff(face, axis=axis) for axis, face in enumerate(result.face_current_a))
        np.testing.assert_allclose(balance, 0, atol=1e-15)
        np.testing.assert_allclose(result.divergence_a_m3, balance / np.prod(spacing), atol=1e-15)
        self.assertLess(abs(sum(result.boundary_current_a.values())), 1e-14)
        self.assertAlmostEqual(result.boundary_power_w / result.dissipated_power_w, 1., places=11)
        self.assertGreater(result.dissipated_power_w, 0)
        self.assertGreaterEqual(result.potential_v.min(), 0)
        self.assertLessEqual(result.potential_v.max(), .02)
        np.testing.assert_array_equal(sigma, original)

    def test_gauge_polarity_and_voltage_controlled_conductivity_scaling(self):
        sigma = np.full((5, 4, 3), .4)
        args = (sigma, (.01, .02, .03))
        base = tank.solve_tank(*args, {"x-": .02, "x+": 0.})
        shifted = tank.solve_tank(*args, {"x-": .07, "x+": .05})
        reverse = tank.solve_tank(*args, {"x-": -.02, "x+": 0.})
        doubled = tank.solve_tank(2 * sigma, args[1], {"x-": .02, "x+": 0.})
        np.testing.assert_allclose(shifted.potential_v, base.potential_v + .05, atol=1e-14)
        np.testing.assert_allclose(shifted.current_density_a_m2, base.current_density_a_m2, atol=1e-12)
        np.testing.assert_allclose(reverse.current_density_a_m2, -base.current_density_a_m2, atol=1e-12)
        np.testing.assert_allclose(doubled.potential_v, base.potential_v, atol=1e-14)
        np.testing.assert_allclose(doubled.current_density_a_m2, 2 * base.current_density_a_m2, atol=1e-12)

    def test_constant_boundary_potential_is_current_off(self):
        result = tank.solve_tank(np.full((3, 4, 2), .4), (.01, .01, .01), {"x-": .02})
        np.testing.assert_allclose(result.potential_v, .02, atol=1e-14)
        np.testing.assert_allclose(result.current_density_a_m2, 0., atol=1e-12)

    def test_nonlinear_harmonic_potential_has_second_order_convergence(self):
        # Independent exact Laplace solution: x^4 - 6*x^2*y^2 + y^4.
        # A linear slab alone cannot measure convergence because it is exact on every mesh.
        errors = []
        for n in (4, 8, 16):
            spacing = (.1 / n,) * 3
            centers = (np.arange(n) + .5) / n
            grid = np.meshgrid(centers, centers, centers, indexing="ij")
            exact = .01 * (grid[0]**4 - 6 * grid[0]**2 * grid[1]**2 + grid[1]**4)
            boundary = {}
            for axis, name in enumerate("xyz"):
                for side, coordinate in (("-", 0.), ("+", 1.)):
                    coords = [centers] * 3
                    coords[axis] = np.array([coordinate])
                    x, y, _ = np.meshgrid(*coords, indexing="ij")
                    boundary[name + side] = np.squeeze(.01 * (x**4 - 6*x*x*y*y + y**4), axis=axis)
            result = tank.solve_tank(np.full((n, n, n), .5), spacing, boundary)
            errors.append(np.linalg.norm(result.potential_v - exact) / np.linalg.norm(exact))
        self.assertGreater(errors[0] / errors[1], 3.4)
        self.assertGreater(errors[1] / errors[2], 3.4)
        self.assertLess(errors[-1], .005)

    def test_invalid_physical_inputs_are_rejected(self):
        sigma = np.ones((3, 3, 3))
        invalid = [
            (sigma * 0, (.01,) * 3, {"x-": 1.}),
            (sigma * np.nan, (.01,) * 3, {"x-": 1.}),
            (sigma.astype(complex) + 1j, (.01,) * 3, {"x-": 1.}),
            (np.ones((3, 3)), (.01,) * 3, {"x-": 1.}),
            (sigma, (0, .01, .01), {"x-": 1.}),
            (sigma, (.01,) * 3, {}),
            (sigma, (.01,) * 3, {"left": 1.}),
            (sigma, (.01,) * 3, {"x-": np.ones((2, 3))}),
            (sigma, (.01,) * 3, {"x-": np.inf}),
        ]
        for args in invalid:
            with self.subTest(args=str(args)[:80]), self.assertRaises(ValueError):
                tank.solve_tank(*args)


if __name__ == "__main__":
    unittest.main()
