# mrs-jax: existing material and integration into MREIT

Audit: **25 October 2026**, [m9h/mrs-jax](https://github.com/m9h/mrs-jax), commit [`3220889281da83f49c618cc5cc623fb6b56b3a33`](https://github.com/m9h/mrs-jax/tree/3220889281da83f49c618cc5cc623fb6b56b3a33). Source was retrieved read-only through GitHub and relevant tests were run from a temporary snapshot. No changes were made to mrs-jax.

The earlier software review omitted this repository. That omission matters: it already contains an MRSI/EPSI forward-model foundation, not just edited single-voxel processing. **Use mrs-jax as our spectroscopy simulation and processing development base**, with FSL-MRS and other established implementations serving as fitting components and independent references where appropriate.

## What is implemented

| Source | Observed implementation | Use in this program |
|---|---|---|
| [`TissueModel` and phantom loader](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/mrsi_sim.py#L45) | Tissue labels, metabolite concentration maps, T1/T2* maps, static B0 offsets and voxel size; FAST segmentation loader with 13 metabolite assignments including `mIns` | Shared anatomical/metabolic phantom for electrical-exposure and MRS-response experiments |
| [`simulate_mrsi`](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/mrsi_sim.py#L207) | Sums concentration-weighted complex basis FIDs, applies T2* decay and static frequency offsets, then spatial FFT | Phase-encoded MRSI forward signal and k-space generation; arbitrary supplied metabolite FIDs can be used |
| [`simulate_mrsi_from_arrays`](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/mrsi_sim.py#L274) | JAX computation through complex signal generation and spatial FFT; one metabolite per call | Gradients through concentrations, T2* and static field parameters; build joint multi-metabolite inference around this path |
| [EPSI functions](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/mrsi_sim.py#L330) | Bipolar trajectory generation, 1D readout simulation, lobe reversal/reconstruction and a JAX array simulator | Existing starting point for time-resolved spectroscopic readout development |
| [`io.py`](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/io.py), [`io_lcmodel.py`](https://github.com/m9h/mrs-jax/blob/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax/io_lcmodel.py) | Siemens TWIX reader via mapVBVD; `MRSData`; LCModel RAW/BASIS readers and writers | Prisma spectroscopy import and exchange with realistic basis sets; validate the actual WIP dimensions and headers |
| [Processing modules](https://github.com/m9h/mrs-jax/tree/3220889281da83f49c618cc5cc623fb6b56b3a33/src/mrs_jax) | MEGA-PRESS, HERMES, JAX spectral registration, preprocessing, phase correction, water-reference calculations and HTML QC | Reuse processing and QC building blocks; preserve the raw complex signal alongside corrected data |

The README roadmap is behind the source: it still lists synthetic MRS, EPSI/MRSI, Philips I/O, LCModel I/O and JAX spectral registration as planned, while relevant implementations exist. Existence of prototype MRSI/EPSI code does not establish a complete whole-brain acquisition/reconstruction pipeline.

## Boundaries that affect MREIT integration

**Time-varying fields are missing.** The MRSI forward model applies `exp(+i 2*pi*delta_f*t)` for a static Hz-valued offset. It does not accept the measured stimulation waveform or a spatially varying `B(r,t)` callback. Add an independently tested phase-history interface evaluated at actual sample times. The analytical MRI examples in this repository use a negative phase convention, so conversion requires an explicit, calibrated convention mapping as well as the Hz/T conversion.

**EPSI is currently simplified.** Spectral evolution is evaluated at the lobe index times the spectral dwell, omitting intra-lobe sample-time evolution. The NumPy version averages concentrations, B0 and T2* over y/z for a 1D projection. It is not a full 3D scanner-timed model. Before fitting small stimulation effects, validate k-space coordinates, centering and reconstruction against a discrete Fourier reference, then include real ADC timing, gradient trajectories, coils and acquisition imperfections.

**The JAX entry points have different capabilities.** `simulate_mrsi_jax` calls NumPy simulation and converts its output to a JAX array. Differentiating through model parameters fails in that wrapper. `simulate_mrsi_from_arrays` is the actual differentiable implementation. The high-level `quantify_mega_press` pipeline also uses NumPy/SciPy operations; the whole package should not be described as end-to-end differentiable without specifying the path.

**Myo-inositol is represented, but its recovery is not established.** The default `make_lorentzian_basis` helper models each named metabolite as one broadened resonance. For quantitative mI response experiments, supply sequence-specific coupled-spin basis FIDs and test concentration recovery under spectral overlap, baseline, linewidth and B0 changes. The current high-level quantification routine targets a GABA Gaussian peak; it is not a general whole-brain mI fitter. The stored T1 map is not consumed by the MRSI signal kernel, so TR saturation and sequence preparation are not modeled there.

## Runtime checks performed

Ran these upstream tests, unmodified, from the pinned snapshot using an existing local environment (NumPy 2.5.4, JAX 0.10.2):

```sh
PYTHONPATH=src python -m pytest \
  tests/test_mrsi_simulator.py \
  tests/test_epsi_readout.py \
  tests/test_simulator_validation.py -q
```

**18 passed, 5 skipped.** The skipped tests require unavailable WAND/ISMRM files at the paths encoded in the tests. This does not verify the README's real-data claims, the entire processing suite, GPU performance or quantitative mI accuracy. Some passing checks establish peak presence or location rather than full reconstruction agreement.

Additional targeted checks found:

- The array-based MRSI concentration gradient agreed with an independent quadratic scaling calculation for a simple constant-basis case.
- Differentiating through `simulate_mrsi_jax` produced `TracerArrayConversionError`, consistent with its NumPy implementation.
- `reconstruct_epsi` changes its input array during in-place lobe reversal. Preserve raw input and remove this side effect before relying on repeatable processing.
- With a basis shorter than the requested EPSI train, NumPy zero-fills later lobes while JAX repeats the last basis sample. Align this behavior and add an equivalence regression before using variable-length basis inputs.

These findings are recorded for subsequent mrs-jax work; its source was not patched in this audit.

## Revised division of responsibilities

| Project | Role |
|---|---|
| **SimNIBS** | Reference anatomy/electrode conductivity solver and E/J/B prediction |
| **sbi4dwi** | Diffusion/microstructure priors, existing inference machinery and proposed scanner/Bloch extensions |
| **mrs-jax** | Metabolic phantoms, basis-FID spatial encoding, MRSI/EPSI forward models and spectroscopy processing |
| **neuroelectric-mri** | Shared physical conventions, cross-simulator benchmarks, Prisma synchronization requirements and reproducible dataset specifications |
| **FSL-MRS / other references** | Realistic basis/fitting components and independent quantitative comparisons, with their own assumptions documented |

The first integration target is a **synthetic disentangling experiment**: vary current-induced phase with concentrations fixed; vary mI concentration with current-induced phase fixed; then vary both. The analysis should distinguish spectral/phase artifacts from concentration changes. Preserve current-off/null cases and raw phase before spectral registration. Parameterize any delayed metabolic response separately from conductivity and current density; do not assume mI is proportional to E or J.

Only after these checks should a shared phantom drive scanner-timed MREIT and MRSI simulations. That reuses the existing code while establishing the additional physics and validation required by the [Prisma program](prisma-platform.md).
