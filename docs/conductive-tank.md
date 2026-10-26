# Conductive-tank forward reference

Implemented **26 October 2026** using [red–green TDD](tdd-conductive-tank.md). This is a small 3D scalar-conductivity reference on an orthogonal Cartesian grid. It solves the steady electrical problem and exports potential, electric field, current density, conservative face currents and power diagnostics. The [synthetic release](../data/synthetic_tank/README.md) includes six controls and analytical convergence measurements.

## Physical and numerical contract

For passive, scalar, positive conductivity, solve `div(sigma grad u) = 0`, with `J = -sigma grad u`. Units are **metres, volts, amperes, S/m, V/m and A/m²**. Arrays use `(x, y, z)` axes, with vector components on the last axis. Cell centers are `(i + 1/2)*h`; the tank spans `[0, nx*hx] × [0, ny*hy] × [0, nz*hz]` in a right-handed frame.

The finite-volume approach enforces flux balance over each cell. [NIST's FiPy numerical documentation](https://pages.nist.gov/fipy/en/master/numerical/discret.html) provides background on cell-centered finite-volume diffusion discretization; this implementation does not use FiPy code and has not been compared numerically with FiPy or SimNIBS.

For adjacent cells p and q, the two half-cell resistances are in series:

```text
R_pq = h/(2*sigma_p*A) + h/(2*sigma_q*A)
G_pq = 1/R_pq
I_pq = G_pq * (u_p - u_q)
```

This harmonic interface conductance is essential for layered materials. A Dirichlet boundary at voltage U contributes a half-cell conductance `G_boundary = 2*sigma*A/h` and outward current `G_boundary*(u_cell-U)`. Each row of the sparse system enforces zero net current from one cell. Unspecified boundary faces are insulated. At least one Dirichlet face fixes the gauge.

The API accepts either a constant potential for an entire face or prescribed values at its face centers. **A varying boundary potential is a mathematical test condition, not a single metal electrode.** The physical plate examples use constant potentials. There are no finite electrode patches, contact impedances, conductivity tensors or zero-conductivity holes in this version.

`solve_plate_current` solves a unit-voltage plate problem and rescales it to the requested total current. Positive current travels from the minus-axis plate to the plus-axis plate, which is grounded. This assumes linear, voltage-independent conductivity and ideal electrodes. It does not prescribe the local injection density: the solution determines that distribution.

## Flux, field and power outputs

- `face_current_a[axis]` is positive in the increasing-coordinate direction on **both** ends of that axis. Its shape has one extra sample in that direction.
- `boundary_current_a` uses outward normals: an injecting electrode has negative total current.
- Cell divergence is the signed sum of its face currents divided by cell volume. These face fluxes, not a differentiated cell-field image, define conservation.
- Cell J averages adjacent face current densities. Cell E is J divided by the cell's scalar conductivity. This is a field reconstruction, not an independent gradient oracle.
- Discrete dissipated power sums `G*(delta_u)**2` over internal and boundary connections, each once. It matches `-sum(U_boundary*I_outward)`. Averaging cell J/E and integrating their product is not the exact discrete power identity.
- Relative boundary imbalance divides net outward current by half the integrated **local absolute normal current**. Taking absolute values only after summing each face can hide simultaneous injection and withdrawal on one face; a regression test explicitly checks this case.

## Independent checks and their limits

The tests verify the homogeneous plate solution in x/y/z, a one-cell half-resistance limit, an aligned layered slab, polarity, gauge shifts, current-off, prescribed total current, per-cell conservation, maximum-principle bounds and discrete power balance. An inclusion exercises heterogeneous conservation; it is not itself independently validated against an anatomical/FEM solver.

The nonlinear convergence oracle is the exact harmonic potential

```text
u = 0.01*(X^4 - 6*X^2*Y^2 + Y^4) volts
X = x/0.1; Y = y/0.1
```

All six boundary faces receive the exact potential. Its continuous Laplacian vanishes, while the discrete solution has measurable truncation error. This makes it a stronger convergence check than a linear slab, which is reproduced to roundoff on every grid. The release includes 4³, 8³, 16³ and 32³ cells; the asymptotic test uses 8³ through 32³. Potential error falls from 0.748% to 0.204% to 0.0526%, with observed orders 1.87 and 1.96. The current-density error is reported separately rather than assumed to have identical convergence.

Under fixed total current, globally multiplying conductivity by three leaves J unchanged and divides potential/E by three. That null is checked on an inhomogeneous tank. Recovering absolute conductivity therefore requires more information than J alone in this ideal-contact model—for example, voltage information or additional justified constraints. It is not legitimate to label an inverse map identifiable solely because it looks plausible.

## Use and next work

```python
import numpy as np
from examples.conductive_tank import solve_plate_current

sigma = np.full((8, 6, 4), 0.5)  # S/m
result = solve_plate_current(sigma, (0.015, 0.08/6, 0.015), 0.001)
print(result.boundary_current_a)  # A, outward sign
print(result.boundary_potential_v["x-"].mean())  # 0.05 V
```

These values define a synthetic example, not equipment settings. A separate [DOLFINx FEM model](fem-complete-electrode.md) now adds finite electrodes, contact impedance and tensor conductivity, including an independent layered-slab comparison with this FV solver. Neither tank exports **B**: an open tank current distribution alone does not specify the complete experimental circuit. Conservative FEM currents, full lead/return geometry, scanner encoding, inversion and acquired-data validation remain later work; see the [FEM roadmap](fem-roadmap.md).
