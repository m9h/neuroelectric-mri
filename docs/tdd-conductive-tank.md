# Red–green TDD record: conductive tank

Executed **26 October 2026**. Tests were written and run in the red state before each implementation. Red and green states are retained as separate commits; only the final passing state is pushed as the branch tip.

| Cycle | Red evidence | Green evidence |
|---|---|---|
| Scalar tank solver | [`9e7d5e3`](https://github.com/m9h/neuroelectric-mri/commit/9e7d5e3): eight test methods fail against a `NotImplementedError` scaffold (18 errors including subtests) | [`b93cc8b`](https://github.com/m9h/neuroelectric-mri/commit/b93cc8b): eight tests pass after sparse finite-volume implementation |
| Fixed-current plates | [`508ad55`](https://github.com/m9h/neuroelectric-mri/commit/508ad55): three new methods fail because `solve_plate_current` is absent; the eight prior methods pass | [`733d980`](https://github.com/m9h/neuroelectric-mri/commit/733d980): all 11 solver methods pass |
| Reproducible release | [`5c330fb`](https://github.com/m9h/neuroelectric-mri/commit/5c330fb): three release methods fail against an unimplemented export function | [`ac079ce`](https://github.com/m9h/neuroelectric-mri/commit/ac079ce): three release methods pass, including repeat generation and hash checks |
| Diagnostic correction | [`2eb1042`](https://github.com/m9h/neuroelectric-mri/commit/2eb1042): a fault-injection ledger expects relative imbalance `0.1/1.05`, but face-net normalization returns `2.0` | [`7a5e862`](https://github.com/m9h/neuroelectric-mri/commit/7a5e862): normalization integrates local absolute boundary currents first; all four release tests pass |

## A test assumption corrected explicitly

The first solver implementation passed seven methods but missed the convergence-rate threshold on the coarsest 4³ grid: its error reduction to 8³ was 3.16, below the requested 3.4. The nonlinear quartic pattern was not yet in the asymptotic range. Refinement gave reductions of 3.67 and 3.88 on 8³→16³→32³, approaching the theoretical factor of four.

The test was changed to use **8³, 16³ and 32³**, retaining its reduction threshold. The 4³ result remains in the published convergence table and plot. No failing numerical result was discarded from the report or replaced with a claim of exact second-order behavior at coarse resolution.

## Reproduction

In a Python environment with `requirements-tank.txt` installed:

```sh
python -m unittest discover -s tests -p 'test_conductive_tank.py' -v
python -m unittest discover -s tests -p 'test_tank_benchmark.py' -v
```

The red commits are deliberately nonworking historical states. Use an isolated checkout to inspect them; do not install them as a usable release. The final workbench contains 37 tests with the complete simulation requirements installed. The new solver/release portion comprises 15 tests; the remaining 22 cover the previous phase, loop and spectroscopy work.

Validation uses analytical potentials/current densities, series resistance, conservation, a discrete power identity, nonlinear mesh refinement and a known conductivity-scale ambiguity. No acquired scanner data, physical equipment, external FEM solver or conductivity inversion was tested in this work package.
