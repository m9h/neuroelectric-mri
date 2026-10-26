import importlib.util
import unittest

AVAILABLE = all(importlib.util.find_spec(p) is not None for p in ("dolfinx", "gmsh", "scipy"))
if AVAILABLE:
    import numpy as np
    import ufl
    from dolfinx import fem
    from examples.fem_tank import build_tank_mesh
    from examples.fem_adaptive import refine_marked, solve_adaptive


@unittest.skipUnless(AVAILABLE, "Optional DOLFINx/Gmsh environment required")
class AdaptiveCEMTests(unittest.TestCase):
    def setUp(self):
        self.currents = {101: .001, 102: -.001}
        self.contacts = {101: .002, 102: .005}

    def test_local_refinement_preserves_material_volume_and_electrode_area(self):
        coarse = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5, layered=True)
        count = coarse.mesh.topology.index_map(3).size_local
        refined = refine_marked(coarse, np.arange(10, dtype=np.int32))
        self.assertGreater(refined.mesh.topology.index_map(3).size_local, count)
        np.testing.assert_array_equal(refined.cell_tags.values,
                                       coarse.cell_tags.values[refined.parent_cell[refined.cell_tags.indices]])
        for data in (coarse, refined):
            dx = ufl.Measure("dx", domain=data.mesh, subdomain_data=data.cell_tags)
            ds = ufl.Measure("ds", domain=data.mesh, subdomain_data=data.facet_tags)
            for tag in (1, 2):
                self.assertAlmostEqual(fem.assemble_scalar(fem.form(1 * dx(tag))), .12 * .08 * .06 / 2, places=13)
            for tag in (101, 102):
                self.assertAlmostEqual(fem.assemble_scalar(fem.form(1 * ds(tag))), .04 * .03, places=13)

    def test_refinement_reduces_energy_gap_and_reports_unmet_target(self):
        data = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5)
        out = solve_adaptive(data, .5, self.currents, self.contacts,
                             relative_gap_tolerance=1e-8, max_refinements=1)
        self.assertFalse(out.converged)
        self.assertEqual(len(out.history), 2)
        before, after = out.history
        self.assertGreater(after["cells"], before["cells"])
        self.assertGreater(after["primal_power_w"], before["primal_power_w"])
        self.assertLess(after["relative_energy_gap"], before["relative_energy_gap"])
        self.assertGreater(after["complementary_power_w"], after["primal_power_w"])
        self.assertAlmostEqual(sum(out.reconstruction.cell_indicator_squared_w), after["estimator_squared_w"], places=12)

    def test_exact_slab_stops_without_refinement(self):
        data = build_tank_mesh()
        out = solve_adaptive(data, .5, self.currents, self.contacts)
        self.assertTrue(out.converged)
        self.assertEqual(len(out.history), 1)
        self.assertIs(out.mesh_data, data)

    def test_invalid_cell_selection_and_adaptive_controls(self):
        data = build_tank_mesh()
        n = data.mesh.topology.index_map(3).size_local
        for cells in ([], [-1], [n], [1.5], [[1]], [True]):
            with self.subTest(cells=cells), self.assertRaises(ValueError):
                refine_marked(data, cells)
        for controls in ({"max_refinements": -1}, {"max_refinements": 1.5},
                         {"relative_gap_tolerance": 0}, {"relative_gap_tolerance": np.nan},
                         {"marking_fraction": 1.1}):
            with self.subTest(controls=controls), self.assertRaises(ValueError):
                solve_adaptive(data, .5, self.currents, self.contacts, **controls)
