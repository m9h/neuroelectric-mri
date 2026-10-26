# 🟢 Adaptive conservative-current benchmark

**Synthetic only.** CG2 complete-electrode potentials, RT degree-2 equilibrated currents and energy-driven mesh refinement. [Numerics and assumptions](../../docs/fem-conservative-current.md) · [TDD record](../../docs/tdd-fem-conservative.md).

![Conservative current, resistance interval and error indicator](summary.png)

The relative complementary-energy gap falls from **2.7893% to 0.8165%** after one refinement. For this 1 mA two-electrode case, the final numerical resistance interval is **73.9996–74.6038 Ω**, equivalent to **73.9996–74.6038 mV**. The 1% energy-gap target is met. This is not a pointwise-current or MRI-phase error bound.

| File | Content |
|---|---|
| `initial.msh` | Original tagged Gmsh/OpenCASCADE tank |
| `mesh.xdmf`, `mesh.h5` | Final refined mesh, material tags and electrode/insulation tags; keep both files together |
| `history.csv` | Cell/dof counts, energy bracket, indicators, current residuals, marking and convergence status per iteration |
| `fields.npz` | Geometry, CG coefficients, RT coefficients/dof maps/orientations, physical current samples, cell indicators and current imbalance |
| `parameters.json` | SI model, numerical controls, basis conventions and explicit stopping status |
| `manifest.json` | Source/artifact hashes, package versions, MPI runtime and platform |

Geometry is in metres, current samples in A/m², cell current imbalance in A, and squared energy indicators in W. The displayed indicator density divides the cell indicator by cell volume; it is an error diagnostic, **not Joule heating**. The old potential-gradient J and new equilibrated J are both sampled at cell centers for comparison.

RT coefficients are basis moments, not point values. They require the stated Basix element, Piola mapping and cell-orientation information. XDMF reload can renumber cells and dofs, so **do not assign the archived coefficient vector directly to an arbitrary reloaded function**. The tested restart route reloads `mesh`, `materials` and `boundaries` from XDMF and re-solves the supplied model; it reproduces the electrode power. Physical cell-center/slice samples can be consumed directly with `allow_pickle=False`.

No leads, B fields, MR observations or conductivity inversion are included. The tests include rotated layered tensors; this released adaptive case uses homogeneous 0.5 S/m conductivity. Repeated numerical generation is tested on a fixed environment; byte identity across platforms is not promised. Generated synthetic numerical artifacts use **CC0-1.0**, with third-party software retaining its own license.
