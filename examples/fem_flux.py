"""Conservative H(div) current reconstruction and CEM energy diagnostics."""


def reconstruct_current(mesh_data, solution, *, degree=1):
    """Reconstruct an equilibrated Raviart–Thomas current in SI units."""
    raise NotImplementedError


def mark_cells(indicator_squared_w, fraction=.5):
    """Select the smallest stable-sorted cell set carrying the given indicator fraction."""
    raise NotImplementedError
