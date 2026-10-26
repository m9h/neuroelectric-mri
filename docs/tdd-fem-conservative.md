# Red–green TDD: conservative currents and adaptive refinement

Executed **26 October 2026**. Each test group was run failing before implementation, with red/green states retained as commits.

| Cycle | Red evidence | Green evidence |
|---|---|---|
| RT reconstruction and indicators | [`b0ffbbc`](https://github.com/m9h/neuroelectric-mri/commit/b0ffbbc): seven methods fail against scaffolds, 13 errors including subtests | [`46d5e62`](https://github.com/m9h/neuroelectric-mri/commit/46d5e62): analytical tensor/slab, conservation, energy identity, polarity, floating electrodes and marking checks pass |
| Adaptive mesh and prescribed-current guard | [`44edce8`](https://github.com/m9h/neuroelectric-mri/commit/44edce8): four adaptive methods fail against scaffolds (14 errors); changed contact data fail to raise the required error | [`a7e2be0`](https://github.com/m9h/neuroelectric-mri/commit/a7e2be0): all eight flux and four adaptive methods pass |
| Adaptive release/restart | [`d3ffae3`](https://github.com/m9h/neuroelectric-mri/commit/d3ffae3): two release methods fail against the exporter scaffold | [`6cadbf9`](https://github.com/m9h/neuroelectric-mri/commit/6cadbf9): hash, field, stopping-status, repeatability and XDMF re-solve checks pass; all 25 FEM methods pass |

Two numerical issues deserve explicit records. First, zero divergence and globally balanced boundary flux do not prove that the requested electrode currents were applied. Doubling stored contact impedance produced a balanced reconstruction carrying half the specified current. The new guard checks each target electrode integral before reconstruction and each reconstructed integral afterward.

Second, DOLFINx's refinement parent array uses pre-reordering child indices. Exposing it directly as a local-cell ancestry array failed for **290 of 2,122 cells** in the layered test. Tag transfer itself was correct; the exported ancestry was wrong. Applying `original_cell_index` fixed it. Material volumes, electrode areas and the ancestry map now pass unchanged assertions.

The first exporter also used an outdated permutation call. Tests exposed the pinned API's distinction between entity quadrature permutations and cell element-DOF permutations; the exporter now records `create_cell_permutations()` / `get_cell_permutation_info()` data. No numerical tolerances were relaxed.

The resulting finite-patch refinement reduces η²/P from 2.7893% to 0.8165%, meets the 1% target, and retains exact-slab and explicitly unmet-target cases. Conservation and energy diagnostics are tested by integrating the returned functions independently, not merely checking their reported metadata.
