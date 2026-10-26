# Synthetic conductive-tank reference

**Numerical forward-model data only.** No scanner measurements, physical current monitoring, conductivity inversion or complete-circuit magnetic field is included. The [solver contract](../../docs/conductive-tank.md) and [red–green TDD record](../../docs/tdd-conductive-tank.md) define what has been checked.

![Potential, current density and mesh convergence](summary.png)

## Reproduce

From the repository root:

```sh
python -m pip install -r requirements-tank.txt
python examples/tank_benchmark.py --output /tmp/tank-reproduction
python -m unittest discover -s tests -p 'test_*tank*.py' -v
```

Add `--plot` to generate the figure when matplotlib is installed. The full `requirements-simulation.txt` includes that dependency through mrs-jax. The tank solver itself does not depend on mrs-jax or JAX.

## Cases and results

The tank is 120 × 80 × 60 mm on an 8 × 6 × 4 grid. Ideal full-face plates drive current along +x, with x+ grounded and the four other faces insulated. Internal faces use harmonic conductances; there is no electrode contact impedance.

| Case | Conductivity | Total current | Required plate voltage |
|---|---|---:|---:|
| Homogeneous | 0.5 S/m | +1 mA | 50.000 mV |
| Layered | 0.2 / 0.8 S/m, equal x lengths | +1 mA | 78.125 mV |
| Inclusion | 0.1 S/m inclusion in 0.5 S/m background | +1 mA | 61.235 mV |
| Scaled inclusion | All inclusion-case conductivities ×3 | +1 mA | 20.412 mV |
| Reverse polarity | Same inclusion | −1 mA | −61.235 mV |
| Current off | Same inclusion | 0 | 0 |

The inclusion occupies indices `[2:6, 1:5, 1:3]`. Its globally scaled version has the same J but one-third of the potential and E. The displayed cyan line is a conductivity contour in the selected slice. The homogeneous and layered voltages also follow independent series-resistance formulas; the inclusion's detailed field does not yet have an independent FEM comparison.

For the separate 100 mm homogeneous cube with nonlinear harmonic boundary data:

| Cells per axis | Relative potential L2 error | Observed order | Relative J L2 error |
|---:|---:|---:|---:|
| 4 | 2.362% | — | 6.190% |
| 8 | 0.748% | 1.659 | 2.090% |
| 16 | 0.204% | 1.875 | 0.641% |
| 32 | 0.0526% | 1.956 | 0.186% |

The convergence example uses mathematical, nonuniform boundary potentials on all faces. It is distinct from the physical equipotential-plate cases. Potential and cell-reconstructed current-density errors are reported separately. Last floating-point digits and plot bytes can vary with NumPy/SciPy/BLAS/matplotlib versions.

## File dictionary

- `fields.npz`: per-case cell conductivity/potential/E/J/divergence and all three face-current arrays, plus SI cell centers/spacing and independent slab potentials. Load with `numpy.load(..., allow_pickle=False)`.
- `metrics.csv`: commanded/recovered total current, plate voltage, boundary and per-cell current imbalance, and discrete power balance for each control.
- `convergence.csv`: analytical relative errors, observed order, current/power diagnostics for every grid.
- `parameters.json`: geometry, coordinate and sign conventions, drive values, boundary conditions and convergence reference.
- `manifest.json`: schema, units via parameters, array shapes/dtypes, dependency versions and SHA-256 hashes of generated files and source.
- `summary.png`: optional static figure, included in the checked-in manifest.

Cell arrays use `(x, y, z)`; E/J have a final 3-component axis. Face currents are **amperes through each face**, positive toward increasing coordinate, including at the tank's minus boundary. `boundary_current_a` in the solver API instead uses outward normals. Cell J/E are reconstructed fields; use face currents for conservation. No B0/scanner transform is implied by this synthetic Cartesian frame.

This generated fixture is dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/); it contains no participant data. Code licensing is separate. A measured release will need its own acquisition provenance and permissions.
