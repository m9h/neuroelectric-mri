"""Energy-indicator-driven conforming mesh refinement for the CEM reference."""

from dataclasses import dataclass
from typing import Any

import numpy as np
from dolfinx import mesh

if __package__:
    from .fem_tank import solve_cem, _scalar
    from .fem_flux import reconstruct_current, mark_cells
else:
    from fem_tank import solve_cem, _scalar
    from fem_flux import reconstruct_current, mark_cells


@dataclass
class RefinedMesh:
    mesh: mesh.Mesh
    cell_tags: mesh.MeshTags
    facet_tags: mesh.MeshTags
    # Maps current local child-cell indices to local cells in the previous mesh.
    parent_cell: np.ndarray


@dataclass
class AdaptiveResult:
    mesh_data: Any
    solution: Any
    reconstruction: Any
    history: list[dict]
    converged: bool


def refine_marked(mesh_data, marked_cells):
    """Refine selected cells while transferring material and electrode tags."""
    domain = mesh_data.mesh
    if domain.comm.size != 1 or domain.topology.cell_type != mesh.CellType.tetrahedron:
        raise ValueError("Local refinement supports serial tetrahedral meshes")
    selected = np.asarray(marked_cells)
    n = domain.topology.index_map(3).size_local
    if (selected.ndim != 1 or not selected.size or selected.dtype.kind not in "iu"
            or np.any(selected < 0) or np.any(selected >= n)):
        raise ValueError("Provide a nonempty one-dimensional integer set of valid cell indices")
    selected = np.unique(selected).astype(np.int32)
    domain.topology.create_entities(1)
    domain.topology.create_connectivity(3, 1)
    edges = mesh.compute_incident_entities(domain.topology, selected, 3, 1)
    refined, parent_cell, parent_facet = mesh.refine(
        domain, edges=edges, option=mesh.RefinementOption.parent_cell_and_facet)
    refined.topology.create_entities(2)
    refined.topology.create_connectivity(2, 3)
    refined.topology.create_connectivity(3, 2)
    cell_tags = mesh.transfer_meshtag(mesh_data.cell_tags, refined, parent_cell)
    facet_tags = mesh.transfer_meshtag(mesh_data.facet_tags, refined, parent_cell, parent_facet)
    # DOLFINx reorders cells when constructing the new mesh. Its transfer
    # routine takes the original refinement arrays; our public ancestry map
    # must instead follow the new mesh's local cell order.
    local_parent = parent_cell[refined.topology.original_cell_index].copy()
    return RefinedMesh(refined, cell_tags, facet_tags, local_parent)


def solve_adaptive(mesh_data, conductivity_s_m, electrode_currents_a,
                   contact_impedance_ohm_m2, *, relative_gap_tolerance=.01,
                   marking_fraction=.5, max_refinements=2, degree=2, flux_degree=2):
    """Refine until eta^2 / supplied_power meets the requested tolerance.

    Each step re-solves CEM and reconstructs an admissible RT current, then
    marks the smallest cell set carrying ``marking_fraction`` of eta^2.
    All edges of selected cells are split, with DOLFINx conforming closure;
    additional cells may therefore refine. This is energy-driven adaptivity,
    not a goal-oriented B/phase estimator. A refinement limit always reports
    ``converged=False`` if the requested tolerance remains unmet.
    """
    tolerance = _scalar(relative_gap_tolerance, "relative_gap_tolerance")
    if tolerance <= 0:
        raise ValueError("relative_gap_tolerance must be positive")
    if (isinstance(max_refinements, (bool, np.bool_))
            or not isinstance(max_refinements, (int, np.integer)) or max_refinements < 0):
        raise ValueError("max_refinements must be a nonnegative integer")
    mark_cells(np.zeros(1), marking_fraction)  # validate even if the first solve is exact
    history = []
    for iteration in range(max_refinements + 1):
        solution = solve_cem(mesh_data, conductivity_s_m, electrode_currents_a,
                              contact_impedance_ohm_m2, degree=degree)
        flux = reconstruct_current(mesh_data, solution, degree=flux_degree)
        power = solution.supplied_power_w
        relative = flux.estimator_squared_w / power if power > 0 else (0. if flux.estimator_squared_w == 0 else float("inf"))
        marked = mark_cells(flux.cell_indicator_squared_w, marking_fraction)
        converged = relative <= tolerance
        history.append(dict(
            iteration=iteration, cells=mesh_data.mesh.topology.index_map(3).size_local,
            potential_dofs=solution.potential_v.function_space.dofmap.index_map.size_local,
            flux_dofs=flux.current_density_a_m2.function_space.dofmap.index_map.size_local,
            primal_power_w=power, complementary_power_w=flux.complementary_power_w,
            estimator_squared_w=flux.estimator_squared_w, relative_energy_gap=relative,
            bulk_mismatch_w=flux.bulk_mismatch_w, contact_mismatch_w=flux.contact_mismatch_w,
            max_cell_current_imbalance_a=float(np.max(abs(flux.cell_current_imbalance_a))),
            relative_constraint_residual=flux.relative_constraint_residual,
            normal_jump_l2_a_m=flux.normal_jump_l2_a_m,
            insulation_leakage_l2_a_m=flux.insulation_leakage_l2_a_m,
            marked_cells=0 if converged or iteration == max_refinements else len(marked),
            converged=converged,
        ))
        if converged or iteration == max_refinements:
            return AdaptiveResult(mesh_data, solution, flux, history, converged)
        mesh_data = refine_marked(mesh_data, marked)
