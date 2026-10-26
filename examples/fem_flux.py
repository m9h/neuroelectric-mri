"""Conservative H(div) current reconstruction and CEM energy diagnostics."""

from dataclasses import dataclass

import numpy as np
import ufl
from dolfinx import fem, mesh
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve

if __package__:
    from .fem_tank import _real, _scalar
else:
    from fem_tank import _real, _scalar


@dataclass
class ConservativeCurrent:
    """Raviart–Thomas field and physical diagnostics in SI units.

    Cell arrays follow DOLFINx cell ordering. The indicator is squared energy
    norm (watts), including the contact contribution on boundary cells. Its
    sum is the primal/complementary power gap for an equilibrated field, up
    to numerical residuals. These are floating-point estimates, not interval
    arithmetic certificates. ``degree`` uses Basix's convention: 1 is the
    lowest RT order (often called RT0 in the literature).
    """

    current_density_a_m2: fem.Function
    degree: int
    cell_current_imbalance_a: np.ndarray
    cell_indicator_squared_w: np.ndarray
    electrode_current_a: dict[int, float]
    divergence_l2_a_m32: float
    normal_jump_l2_a_m: float
    insulation_leakage_l2_a_m: float
    bulk_mismatch_w: float
    contact_mismatch_w: float
    estimator_squared_w: float
    complementary_power_w: float
    relative_constraint_residual: float
    relative_stationarity_residual: float


def _matrix(form):
    assembled = fem.assemble_matrix(fem.form(form))
    assembled.scatter_reverse()
    # The copy owns its storage after the DOLFINx assembly object is released.
    return assembled.to_scipy().copy()


def _integral(form):
    return float(fem.assemble_scalar(fem.form(form)))


def _cell_integrals(domain, form_factory):
    space = fem.functionspace(domain, ("DG", 0))
    test = ufl.TestFunction(space)
    vector = fem.assemble_vector(fem.form(form_factory(test)))
    n = domain.topology.index_map(3).size_local
    return vector.array[np.array([space.dofmap.cell_dofs(c)[0] for c in range(n)])].copy()


def reconstruct_current(mesh_data, solution, *, degree=1):
    """Minimize the conductivity-weighted current correction under equilibrium.

    Use RT(degree)/DG(degree-1) on serial affine tetrahedra. Project the CEM
    outward contact flux (u-U)/z onto the RT boundary normal-trace space;
    insulation has exactly zero normal trace. Preserve each electrode's
    total current, then enforce div(J)=0 in the complete divergence space.

    Among fields with that boundary trace and divergence, minimize
    integral (J + sigma grad(u)) . inv(sigma) (J + sigma grad(u)). A global
    saddle-point solve is the reference implementation, not a scalable local
    patch algorithm. One redundant DG multiplier equation is removed; the
    full divergence residual (including that equation) is checked afterward.

    Reconstructed current must not be interpreted as an independently solved
    constitutive field: its difference from -sigma grad(u), including the
    contact-law difference, is the error indicator being reported.
    """
    domain = mesh_data.mesh
    if (domain is not solution.potential_v.function_space.mesh
            or domain is not solution.conductivity_s_m.function_space.mesh):
        raise ValueError("Solution and tags must belong to the same mesh object")
    if domain.comm.size != 1 or domain.topology.dim != 3 or domain.geometry.dim != 3:
        raise ValueError("This reconstruction requires a serial three-dimensional mesh")
    if domain.topology.cell_type != mesh.CellType.tetrahedron or domain.geometry.cmaps[0].degree != 1:
        raise ValueError("Only affine tetrahedra are supported")
    if isinstance(degree, (bool, np.bool_)) or not isinstance(degree, (int, np.integer)) or degree not in (1, 2):
        raise ValueError("RT degree must be integer 1 or 2")
    if solution.relative_residual > 1e-8:
        raise ValueError("The input CEM solution is not sufficiently converged")
    contacts = solution.contact_impedance_ohm_m2
    tags = sorted(solution.applied_current_a)
    if set(contacts) != set(tags) or any(not np.isfinite(z) or z <= 0 for z in contacts.values()):
        raise ValueError("Positive finite contact impedances are required for every electrode")
    domain.topology.create_connectivity(2, 3)
    exterior = mesh.exterior_facet_indices(domain.topology)
    active_facets = np.concatenate([mesh_data.facet_tags.find(t) for t in tags])
    insulation_tags = mesh.meshtags(domain, 2, exterior,
        np.isin(exterior, active_facets).astype(np.int32))
    dx = ufl.Measure("dx", domain=domain)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=mesh_data.facet_tags)
    insulated = ufl.Measure("ds", domain=domain, subdomain_data=insulation_tags)(0)
    n = ufl.FacetNormal(domain)
    potential, sigma = solution.potential_v, solution.conductivity_s_m
    inverse_sigma = ufl.inv(sigma)
    flux_space = fem.functionspace(domain, ("RT", degree))
    multiplier_space = fem.functionspace(domain, ("DG", degree - 1))
    trial, test = ufl.TrialFunction(flux_space), ufl.TestFunction(flux_space)
    eta = ufl.TestFunction(multiplier_space)

    # A normal-trace mass projection handles facet orientations and all
    # moments of each RT order. Boundary constants are reproduced exactly.
    trace_mass = _matrix(ufl.dot(trial, n) * ufl.dot(test, n) * ds)
    trace_rhs_form = sum((potential - solution.electrode_potential_v[t]) / contacts[t]
                         * ufl.dot(test, n) * ds(t) for t in tags)
    trace_rhs = fem.assemble_vector(fem.form(trace_rhs_form)).array.copy()
    boundary_dofs = fem.locate_dofs_topological(flux_space, 2, exterior)
    size = flux_space.dofmap.index_map.size_local
    coefficients = np.zeros(size)
    coefficients[boundary_dofs] = spsolve(trace_mass[boundary_dofs][:, boundary_dofs], trace_rhs[boundary_dofs])
    free = np.setdiff1d(np.arange(size), boundary_dofs)
    mass = _matrix(ufl.inner(ufl.dot(inverse_sigma, trial), test) * dx)
    divergence = _matrix(ufl.div(trial) * eta * dx)
    rhs = fem.assemble_vector(fem.form(-ufl.inner(ufl.grad(potential), test) * dx)).array.copy()
    # DG Lagrange bases form a partition of unity; one multiplier coefficient
    # fixes its constant nullspace. Its omitted constraint is not discarded
    # from the physical residual checks below.
    pressure_free = np.arange(multiplier_space.dofmap.index_map.size_local - 1)
    coupling = divergence[pressure_free][:, free]
    matrix = bmat([[mass[free][:, free], coupling.T], [coupling, None]], format="csr")
    constrained_rhs = np.concatenate([
        rhs[free] - mass[free][:, boundary_dofs] @ coefficients[boundary_dofs],
        -divergence[pressure_free][:, boundary_dofs] @ coefficients[boundary_dofs],
    ])
    result = spsolve(matrix, constrained_rhs)
    if not np.all(np.isfinite(result)) or not np.all(np.isfinite(coefficients)):
        raise RuntimeError("Conservative reconstruction produced nonfinite values")
    coefficients[free] = result[:len(free)]
    current_scale = .5 * sum(abs(i) for i in solution.applied_current_a.values())
    constraint_residual = float(np.linalg.norm(divergence @ coefficients) / (current_scale or 1.))
    stationarity_residual = float(np.linalg.norm(matrix @ result - constrained_rhs)
                                  / (np.linalg.norm(constrained_rhs) or 1.))
    if max(constraint_residual, stationarity_residual) > 1e-8:
        raise RuntimeError("Conservative reconstruction failed its full equilibrium residual check")
    current = fem.Function(flux_space, name="J_equilibrated_A_per_m2")
    current.x.array[:] = coefficients
    current.x.scatter_forward()

    mismatch = current + ufl.dot(sigma, ufl.grad(potential))
    bulk_density = ufl.inner(ufl.dot(inverse_sigma, mismatch), mismatch)
    normal = ufl.dot(current, n)
    contact_density = {t: contacts[t] * (normal - (potential - solution.electrode_potential_v[t])
                                                    / contacts[t])**2 for t in tags}
    bulk_error = _integral(bulk_density * dx)
    contact_error = sum(_integral(contact_density[t] * ds(t)) for t in tags)
    dual_power = _integral(ufl.inner(ufl.dot(inverse_sigma, current), current) * dx)
    dual_power += sum(_integral(contacts[t] * normal**2 * ds(t)) for t in tags)
    indicator = _cell_integrals(domain, lambda q: bulk_density * q * dx
                                + sum(contact_density[t] * q * ds(t) for t in tags))
    cell_balance = _cell_integrals(domain, lambda q: ufl.div(current) * q * dx)
    electrode_current = {t: -_integral(normal * ds(t)) for t in tags}
    return ConservativeCurrent(
        current, degree, cell_balance, indicator, electrode_current,
        np.sqrt(max(0., _integral(ufl.div(current)**2 * dx))),
        np.sqrt(max(0., _integral(ufl.jump(current, n)**2 * ufl.dS))),
        np.sqrt(max(0., _integral(normal**2 * insulated))),
        bulk_error, contact_error, bulk_error + contact_error, dual_power,
        constraint_residual, stationarity_residual,
    )


def mark_cells(indicator_squared_w, fraction=.5):
    """Select the smallest stable-sorted cell set carrying the given indicator fraction."""
    values = _real(indicator_squared_w, "indicator_squared_w")
    fraction = _scalar(fraction, "fraction")
    if values.ndim != 1 or np.any(values < 0) or not 0 < fraction <= 1:
        raise ValueError("Need a nonnegative 1D indicator and 0 < fraction <= 1")
    total = float(np.sum(values))
    if not np.isfinite(total):
        raise ValueError("Indicator sum must be finite")
    order = np.argsort(-values, kind="stable")
    order = order[values[order] > 0]
    if total == 0 or fraction == 1:
        return order
    count = np.searchsorted(np.cumsum(values[order]), fraction * total, side="left") + 1
    return order[:count]
