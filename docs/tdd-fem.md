# Red–green TDD: complete-electrode FEM

Executed **26 October 2026**. Tests were run failing before implementation; the commits preserve that history.

| Cycle | Red | Green |
|---|---|---|
| DOLFINx CEM and tagged CAD | [`37f1345`](https://github.com/m9h/neuroelectric-mri/commit/37f1345): seven methods fail against `NotImplementedError` scaffolds (12 errors including subtests) | [`12c3b94`](https://github.com/m9h/neuroelectric-mri/commit/12c3b94): all seven pass, including analytical/FV comparisons and rotated tensors |
| FEM release and refinement | [`1b5982f`](https://github.com/m9h/neuroelectric-mri/commit/1b5982f): three release methods fail against an unimplemented exporter (seven errors including subtests) | [`70746c6`](https://github.com/m9h/neuroelectric-mri/commit/70746c6): all ten FEM methods pass with a working exporter and refined electrode edges |
| CAD electrode identity | [`5f18cf2`](https://github.com/m9h/neuroelectric-mri/commit/5f18cf2): an electrode covering half a face is incorrectly tagged as the whole face, giving 0.0048 instead of 0.0024 m² | [`da5f81f`](https://github.com/m9h/neuroelectric-mri/commit/da5f81f): preserve electrode ancestry through OpenCASCADE Boolean fragments; all eleven FEM methods pass |

## A convergence failure drove a mesh improvement

The first working exporter passed repeatability and input validation, but failed the prewritten **1% successive CG2 voltage-change criterion**. Uniform target mesh sizes 30/20/12 mm gave:

| Target h | Tetrahedra | CG2 voltage |
|---|---:|---:|
| 30 mm | 194 | 70.2689 mV |
| 20 mm | 534 | 70.9495 mV |
| 12 mm | 1,922 | 72.7641 mV |

The last change was **2.4938%**. Introducing a distance-based mesh field with h/4 at electrode edges reduced it to **1.0384%**, still failing. Refining the edges to h/8 gave **0.3470%**, with 1,238/1,937/4,462 tetrahedra. The threshold and original three target sizes were retained. The initial numerical results are recorded above; the released meshes use the final h/8 grading.

This checks stabilization of electrode voltage. It does not certify local current accuracy or prove a formal asymptotic convergence rate near contact/insulation transitions. CG1/CG2 are nested on each individual mesh; the independently generated h levels are not nested meshes.

## A CAD ambiguity caught during review

With patch fraction √0.5, the electrode and its surrounding insulating face have identical areas and centroids. The initial geometric labeling heuristic therefore selected both. A new regression test demonstrated the doubled area. The fix uses the original electrode's OpenCASCADE Boolean-fragment mapping, preserving topological identity rather than guessing from summary geometry. The default quarter-area release is unchanged numerically.

Local validation: eleven FEM tests pass in the DOLFINx environment; the previous 37 tests pass in the MRS/JAX environment, where optional FEM tests skip. The native Gmsh library initially reported sandbox-denied host-name calls; the successful FEM runs were repeated with those system calls permitted. No test failure was resolved by loosening numerical acceptance thresholds.
