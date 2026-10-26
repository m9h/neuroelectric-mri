"""Energy-indicator-driven conforming mesh refinement for the CEM reference."""


def refine_marked(mesh_data, marked_cells):
    """Refine selected cells while transferring material and electrode tags."""
    raise NotImplementedError


def solve_adaptive(mesh_data, conductivity_s_m, electrode_currents_a,
                   contact_impedance_ohm_m2, *, relative_gap_tolerance=.01,
                   marking_fraction=.5, max_refinements=2, degree=2, flux_degree=2):
    """Refine until the relative complementary-energy gap meets the target or limit."""
    raise NotImplementedError
