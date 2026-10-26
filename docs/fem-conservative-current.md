# Conservative current and adaptive CEM error control

Implemented **26 October 2026** in [fem_flux.py](../examples/fem_flux.py), [fem_adaptive.py](../examples/fem_adaptive.py) and the [adaptive release](../data/synthetic_fem_adaptive/README.md). The original CG potential-gradient current remains available for comparison; the new current is H(div)-conforming and equilibrated.

## Reconstruction

For the solved CEM potential uₕ and electrode voltages Uₕ, find J* in a Raviart–Thomas space that minimizes

\[
\int_\Omega (J^*+\sigma\nabla u_h)^T\sigma^{-1}(J^*+\sigma\nabla u_h)\,dV
\]

subject to zero divergence and the prescribed boundary normal trace. On electrodes, project (uₕ−Uₕ)/z into the RT normal-trace space; on insulation set it to zero. That projection preserves the integral on every electrode. A boundary trace-mass solve handles orientation and all facet moments rather than assuming a particular RT coefficient sign.

The implementation supports **Basix RT degrees 1 and 2**, paired with DG0/DG1 divergence multipliers. Basix degree 1 is the lowest RT order, often called RT0 in the literature. See the [official DOLFINx mixed-Poisson example](https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_mixed-poisson.html). The global saddle-point solve uses SciPy in serial on affine tetrahedra. One multiplier coefficient fixes its constant gauge; the full divergence residual, including the omitted redundant equation, is checked afterward.

Tests integrate divergence in every cell, evaluate normal jumps across interior faces, verify insulation and electrode currents, and recover known homogeneous/layered currents with rotated conductivity tensors. A fault test changes contact parameters so that a globally balanced trace carries only half the specified current. Reconstruction now rejects this: global balance alone is insufficient.

## Energy identity and its scope

Let Jₕ=−σ∇uₕ and gₕ=(uₕ−Uₕ)/z. The squared indicator is

\[
\eta^2=\int_\Omega (J^*-J_h)^T\sigma^{-1}(J^*-J_h)\,dV
+\sum_e z_e\int_{E_e}(J^*\cdot n-g_h)^2\,dS.
\]

The contact term matters: projecting a quadratic boundary potential onto a lower-degree normal trace loses information, even while preserving total current. It is included in every indicator and assigned to the adjacent boundary cell.

For an admissible current, define complementary power

\[
D(J^*)=\int_\Omega J^{*T}\sigma^{-1}J^*\,dV
+\sum_e z_e\int_{E_e}(J^*\cdot n)^2\,dS.
\]

Integration by parts, zero divergence and the prescribed electrode totals give **η² = D(J*) − Pₕ**, where Pₕ=ΣₑUₕ,ₑIₑ for the converged Galerkin solution. The exact CEM power satisfies **Pₕ ≤ P ≤ D(J*)**. Equivalently, η² is the sum of the squared primal energy error and the squared dual current/contact error. This is the complementary-energy argument specialized here to our inward-positive current convention and SPD conductivity. [Díaz-Avalos (2023)](https://doi.org/10.17268/sel.mat.2023.01.09) develops CEM duality and a posteriori estimates; its publisher abstract/metadata were inspected as background, not its numerical implementation reproduced.

For a two-electrode current ±I, divide the powers by I² to obtain a resistance interval. The code checks the power identity, full equilibrium residual and individual electrode currents before reporting the estimate. Results are **floating-point numerical brackets under these assumptions**, not interval-arithmetic certificates. They bound this CEM energy/resistance problem, not geometry error, uncertain material/contact values, pointwise J, or MRI phase. A 1% relative energy gap also does not mean 1% relative energy-norm error: the norm involves a square root.

## Adaptive loop

The loop solves CG2 CEM, reconstructs RT degree 2, and marks the smallest stable-sorted set of cells carrying 50% of η². All edges of those cells are selected for DOLFINx refinement, with conforming closure. Cell/facet tags are transferred through the parent maps. DOLFINx reorders child cells; the public ancestry array is explicitly permuted into the new local cell ordering.

Stop when **η²/Pₕ ≤ 0.01** or the configured refinement limit is reached. An unmet target returns `converged=False` and `stop_reason=refinement_limit` in the release. The zero-error slab stops without refining. This is implemented energy-driven adaptivity; B/phase goal-oriented adaptivity remains future work.

The released finite-patch case gives:

| Level | Tetrahedra | CG2 dofs | RT dofs | Relative energy gap | Resistance interval |
|---|---:|---:|---:|---:|---:|
| Initial | 1,238 | 2,381 | 12,168 | 2.7893% | 73.4761–75.5255 Ω |
| Adaptive | 6,826 | 10,926 | 63,921 | **0.8165%** | **73.9996–74.6038 Ω** |

Only 93 initial cells were marked, but conforming closure increased the total to 6,826. This overhead is recorded rather than hidden. The final maximum cell current imbalance is about **4.05×10⁻¹⁷ A** at 1 mA applied current. These results do not establish scaling to anatomical meshes; local patch reconstruction, better refinement efficiency and PETSc preconditioning remain opportunities.

## Reproduce and next gate

In the [isolated FEM environment](fem-complete-electrode.md#reproduce):

```sh
python -m unittest discover -s tests -p 'test_fem*.py' -v
python examples/fem_adaptive_benchmark.py --output /tmp/synthetic-fem-adaptive --plot
```

There are 25 FEM tests, including 14 added for reconstruction, refinement and release/restart behavior. With the 37 spectroscopy/FV tests, this is **62 distinct tests** across the two environments. [TDD history](tdd-fem-conservative.md).

The next external validation gate is an independent GetDP CEM solve on the same tagged mesh, followed by outgoing/return-lead geometry and complete-circuit B. Broader manufactured tensor fields, mesh-quality/geometry-order studies, B-oriented estimators, MPI scaling and adjoints remain open in the [FEM roadmap](fem-roadmap.md).
