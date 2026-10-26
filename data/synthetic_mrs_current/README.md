# Synthetic current-phase / concentration controls

**Software benchmark only. No acquired MRI, measured current, conductivity inference or quantitative mI validation.** The `mIns` basis is a single illustrative Lorentzian, and all fitting uses known T2*/B0 and the same basis as simulation. "Oracle" fitting also knows the injected phase.

The phantom contains 6 × 6 × 1 voxels and 256 FID samples. A complete square-loop Biot–Savart field drives a piecewise current waveform, which mrs-jax integrates at each sample before spatial encoding. Each phase encode assumes an identical excitation-relative waveform. The central 2 × 2 region receives a separate nominal mI concentration change in the relevant controls.

## Reproduce

From the repository root, in an environment with the pinned simulation dependency:

```sh
python -m pip install -r requirements-simulation.txt
python examples/mrs_current_benchmark.py --output /tmp/mrs-current-reproduction
python -m unittest discover -s tests -v
```

Optional synthetic k-space noise is enabled with `--noise-sd 0.01 --seed 19`. Noise SD applies separately to each real and imaginary k-space component; it is not a calibrated scanner noise level. The checked-in fixture is noiseless. Floating-point last digits can depend on NumPy/BLAS versions; the tests use numerical tolerances.

## Results

ROI change in **nominal mM under the illustrative basis**, relative to 7:

| Control | Injected change | Fitting without current phase | Fitting with known phase |
|---|---:|---:|---:|
| Null | 0 | 0.0000 | 0.0000 |
| Phase only | 0 | −0.0417 | 0.0000 |
| Concentration only | +0.3500 | +0.3500 | +0.3500 |
| Combined | +0.3500 | +0.3060 | +0.3500 |
| Reverse polarity | 0 | −0.0436 | 0.0000 |
| Current delayed by 3 ms | 0 | −0.0287 | 0.0000 |

In the delayed case, correcting with the wrong zero-delay waveform leaves a +0.0066 nominal mM apparent change. These numbers illustrate model mismatch in this one synthetic setting; they do not predict a human stimulation effect. Matched-model noiseless recovery is a necessary implementation check and is not independent scientific validation.

## Files and conventions

- `signals.npz`: complex k-space per control; injected voxel phase and concentrations; fitted maps; field per ampere, geometry coordinates, basis FIDs, sample times and ROI. Load with `numpy.load(..., allow_pickle=False)`.
- `parameters.json`: SI current/timing/geometry specification, basis settings, metabolite order, sign convention and seed.
- `metrics.csv`: unrounded results, including the incorrect zero-delay fit.
- `manifest.json`: array shapes/dtypes, file and source SHA-256 hashes, mrs-jax revision/worktree state and Python/NumPy/SciPy versions.

`*_kspace` has axes `(kx, ky, kz, time)` in NumPy FFT ordering. `*_phase_rad` has `(x, y, z, time)`; concentration and fit arrays have `(x, y, z, metabolite)`. The final axis follows `metabolite_order`. `field_per_amp_t` is B0-parallel T/A; `positions_m` uses a right-handed synthetic scanner frame, B0 along +z. The signal's additional phase is `+gamma integral I*q*B_per_amp dt`; negate the separate analytical workbench's negative-sign phase when comparing conventions. Chemical-shift/static-B0 evolution remains in the MRS kernel.

This small generated fixture is dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). It contains no participant data. Reuse of the generator and dependency source follows each repository's separate code licensing terms. A measured release will need its own dataset card and permissions.
