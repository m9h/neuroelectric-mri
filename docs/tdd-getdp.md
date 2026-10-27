# Red–green TDD: independent GetDP CEM

Executed **26 October 2026**. The failing tests were run and committed before implementing each stage.

| Cycle | Red evidence | Green evidence |
|---|---|---|
| Native independent CEM | [`83f0e83`](https://github.com/m9h/neuroelectric-mri/commit/83f0e83): six test methods fail against the scaffold, 16 errors including subtests | [`0ff5cce`](https://github.com/m9h/neuroelectric-mri/commit/0ff5cce): analytical P1/P2 slabs, layered resistance, tensor conductivity, contact power, reciprocity, gauge, polarity, zero current, invalid inputs and retained native outputs pass |
| Independent comparison and release | [`41d879f`](https://github.com/m9h/neuroelectric-mri/commit/41d879f): three test methods fail, six errors including subtests | [`1a68f8d`](https://github.com/m9h/neuroelectric-mri/commit/1a68f8d): same-mesh P1/P2 scalar/tensor fields, geometry/tag/current-sign round-trip, nine-case export, artifact hashes and adaptive energy interval pass |

The first native run exposed two GetDP integration details. `Order` is reserved syntax, so the generated degree variable is named `TissueDegree`. More substantially, an integral postprocessing quantity must specify its integration group when printed globally: `Print[quantity[Region], OnGlobal, ...]`. Omitting the group returned zeros for voltage and power while local potential samples were nonzero. The analytical voltage and independently integrated power/current tests caught this immediately. No assertions or tolerances were relaxed.

The mesh exporter checks *geometric* connectivity and labels after import, rather than comparing local node IDs. It uses DOLFINx's topology-to-geometry map and positive tetrahedral orientation; the native GetDP re-solve then checks current direction and the electrode voltage. This is necessary because XDMF/DOLFINx reordering need not preserve entity IDs.

Locally, all nine new tests pass with GetDP 4.0.0-git-259f4c38. All 25 prior FEM tests and all 37 spectroscopy/FV tests also pass in their respective environments. Linux CI [run 40455074002](https://github.com/m9h/neuroelectric-mri/actions/runs/40455074002) passes the same physics/comparison assertions with the checksum-pinned GetDP 3.5.0 binary. The synthetic job discovers 71 tests and skips the 34 native FEM/GetDP methods; the FEM job executes those 34 without skips. This gives **71 distinct passing tests**, not 105 independent tests.
