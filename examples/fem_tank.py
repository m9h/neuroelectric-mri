"""DOLFINx complete-electrode tank benchmark (SI units)."""


def build_tank_mesh(lengths_m=(.12, .08, .06), mesh_size_m=.02,
                    patch_fraction=1., layered=False, mesh_path=None):
    """Create an OpenCASCADE tank with tagged materials and electrodes."""
    raise NotImplementedError


def solve_cem(mesh_data, conductivity_s_m, electrode_currents_a,
              contact_impedance_ohm_m2, *, degree=1, reference_tag=None):
    """Solve the current-driven complete electrode model."""
    raise NotImplementedError
