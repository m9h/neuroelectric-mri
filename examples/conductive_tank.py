"""Conservative 3D scalar-conductivity tank reference, in SI units.

Solve div(sigma grad u)=0 on an orthogonal cell-centered Cartesian grid.
Internal faces use two half-cell resistances in series (harmonic conductivity).
Specified boundary faces have Dirichlet potentials; other faces are insulated.
This is a small forward reference, not an anatomical solver or MREIT inversion.

Finite-volume background: https://pages.nist.gov/fipy/en/master/numerical/discret.html
The implementation here is independent of FiPy; analytical tests are the oracle.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


@dataclass
class TankSolution:
    """Potential, conservative fluxes and diagnostic fields for a tank.

    All cell arrays have axes (x, y, z); vector components are the final axis.
    face_current_a[axis] has one extra entry along that axis, with positive
    current toward +axis, including at both exterior faces. Boundary totals
    instead use the outward-normal convention: injection is negative.

    Cell J is the arithmetic average of adjacent face current densities.
    Cell E is J/sigma. These reconstructions do not replace the conservative
    face currents for divergence or the face-resistance power calculation.
    """

    potential_v: np.ndarray
    conductivity_s_m: np.ndarray
    spacing_m: np.ndarray
    boundary_potential_v: dict[str, np.ndarray]
    face_current_a: tuple[np.ndarray, np.ndarray, np.ndarray]
    current_density_a_m2: np.ndarray
    electric_field_v_m: np.ndarray
    divergence_a_m3: np.ndarray
    boundary_current_a: dict[str, float]
    boundary_power_w: float
    dissipated_power_w: float


def _real_array(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _inputs(conductivity_s_m, spacing_m, boundary_potential_v):
    sigma = _real_array(conductivity_s_m, "conductivity_s_m")
    spacing = _real_array(spacing_m, "spacing_m")
    if sigma.ndim != 3 or not sigma.size or np.any(sigma <= 0):
        raise ValueError("conductivity_s_m must be a nonempty 3D positive scalar map")
    if spacing.shape != (3,) or np.any(spacing <= 0):
        raise ValueError("spacing_m must contain three positive cell sizes")
    if not isinstance(boundary_potential_v, Mapping) or not boundary_potential_v:
        raise ValueError("At least one Dirichlet face is required to fix the potential gauge")
    boundaries = {}
    for name, value in boundary_potential_v.items():
        if name not in ("x-", "x+", "y-", "y+", "z-", "z+"):
            raise ValueError(f"Unknown boundary face: {name}")
        axis = "xyz".index(name[0])
        shape = sigma.shape[:axis] + sigma.shape[axis + 1:]
        value = _real_array(value, f"Boundary {name}")
        if value.shape not in ((), shape):
            raise ValueError(f"Boundary {name} needs a scalar or transverse shape {shape}")
        boundaries[name] = np.broadcast_to(value, shape).copy()
    return sigma, spacing, boundaries


def _slices(axis):
    low, high = [slice(None)] * 3, [slice(None)] * 3
    low[axis], high[axis] = slice(None, -1), slice(1, None)
    return tuple(low), tuple(high)


def _conductances(sigma, spacing):
    conductances = []
    for axis, name in enumerate("xyz"):
        area = np.prod(np.delete(spacing, axis))
        shape = list(sigma.shape)
        shape[axis] += 1
        g = np.empty(shape)
        interior = [slice(None)] * 3
        interior[axis] = slice(1, -1)
        low, high = _slices(axis)
        # A / (h/(2*sigma_left) + h/(2*sigma_right)), in siemens.
        g[tuple(interior)] = area / (spacing[axis] / (2 * sigma[low])
                                     + spacing[axis] / (2 * sigma[high]))
        for side, index in (("-", 0), ("+", -1)):
            face = [slice(None)] * 3
            face[axis] = index
            g[tuple(face)] = 2 * sigma[tuple(face)] * area / spacing[axis]
        if not np.all(np.isfinite(g)) or np.any(g <= 0):
            raise ValueError("Conductance is not finite and positive; rescale the inputs")
        conductances.append(g)
    return tuple(conductances)


def _recover(potential, sigma, spacing, boundaries, conductances):
    currents, density = [], []
    boundary_current = {}
    power_boundary, power_dissipated = 0., 0.
    for axis, name in enumerate("xyz"):
        g = conductances[axis]
        current = np.zeros_like(g)
        low, high = _slices(axis)
        interior = [slice(None)] * 3
        interior[axis] = slice(1, -1)
        interior = tuple(interior)
        du = potential[low] - potential[high]
        current[interior] = g[interior] * du
        power_dissipated += float(np.sum(g[interior] * du**2))
        for side, index, outward_sign in (("-", 0, -1), ("+", -1, 1)):
            key = name + side
            face = [slice(None)] * 3
            face[axis] = index
            face = tuple(face)
            boundary_current[key] = 0.
            if key in boundaries:
                du_out = potential[face] - boundaries[key]
                current_out = g[face] * du_out
                current[face] = outward_sign * current_out
                boundary_current[key] = float(np.sum(current_out))
                power_boundary -= float(np.sum(boundaries[key] * current_out))
                power_dissipated += float(np.sum(g[face] * du_out**2))
        currents.append(current)
        area = np.prod(np.delete(spacing, axis))
        density.append(.5 * (current[low] + current[high]) / area)
    j = np.stack(density, axis=-1)
    divergence = sum(np.diff(f, axis=a) for a, f in enumerate(currents)) / np.prod(spacing)
    return TankSolution(potential, sigma, spacing, boundaries, tuple(currents),
                        j, j / sigma[..., None], divergence, boundary_current,
                        power_boundary, power_dissipated)


def solve_tank(
    conductivity_s_m: np.ndarray,
    spacing_m: tuple[float, float, float],
    boundary_potential_v: Mapping[str, float | np.ndarray],
) -> TankSolution:
    """Solve a passive rectangular tank using a sparse finite-volume system.

    Parameters
    ----------
    conductivity_s_m : array, shape (nx, ny, nz)
        Positive, finite scalar conductivity at cell centers, in S/m.
        Tensor conductivity, zero-conductivity holes and contact impedance
        are not supported. Interfaces lie on cell faces.
    spacing_m : sequence of three floats
        Uniform cell size along each axis, in metres. Domain is [0, n*h].
    boundary_potential_v : mapping
        Keys x-/x+/y-/y+/z-/z+ set entire face potentials in volts. Values
        are scalars or arrays with the remaining two axes in their original
        order, sampled at face centers. Unspecified faces have zero normal
        current. A nonuniform face value is mathematical boundary data, not
        a single equipotential metal electrode.

    Returns
    -------
    TankSolution
        Potential, face currents, reconstructed E/J and conservation/power
        diagnostics. No magnetic field, leads or conductivity inversion.
    """
    sigma, spacing, boundaries = _inputs(conductivity_s_m, spacing_m, boundary_potential_v)
    conductances = _conductances(sigma, spacing)
    ids = np.arange(sigma.size).reshape(sigma.shape)
    diagonal = np.zeros(sigma.size)
    rhs = np.zeros(sigma.size)
    rows, columns, values = [], [], []
    for axis, name in enumerate("xyz"):
        low, high = _slices(axis)
        p, q = ids[low].ravel(), ids[high].ravel()
        interior = [slice(None)] * 3
        interior[axis] = slice(1, -1)
        g = conductances[axis][tuple(interior)].ravel()
        diagonal[p] += g
        diagonal[q] += g
        rows.extend([p, q])
        columns.extend([q, p])
        values.extend([-g, -g])
        for side, index in (("-", 0), ("+", -1)):
            key = name + side
            if key in boundaries:
                face = [slice(None)] * 3
                face[axis] = index
                face = tuple(face)
                p = ids[face].ravel()
                g = conductances[axis][face].ravel()
                diagonal[p] += g
                rhs[p] += g * boundaries[key].ravel()
    rows.append(ids.ravel())
    columns.append(ids.ravel())
    values.append(diagonal)
    matrix = coo_matrix((np.concatenate(values),
                         (np.concatenate(rows), np.concatenate(columns))),
                        shape=(sigma.size, sigma.size)).tocsr()
    potential = spsolve(matrix, rhs).reshape(sigma.shape)
    if not np.all(np.isfinite(potential)):
        raise RuntimeError("Tank solve produced nonfinite potentials")
    return _recover(potential, sigma, spacing, boundaries, conductances)
