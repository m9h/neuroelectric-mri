"""Independent GetDP complete-electrode reference (SI units)."""


def solve_getdp(mesh_path, conductivity_s_m, electrode_currents_a,
                contact_impedance_ohm_m2, *, degree=1, reference_tag=None,
                output_dir=None, executable=None):
    """Run GetDP's own CEM assembly/solve and return its fields and diagnostics."""
    raise NotImplementedError
