# Implementation roadmap: models → instruments → open measurements

Updated **26 October 2026**. Scanner target: **Siemens Prisma 3 T**, with WIPs and a route to Siemens development. Progress is gated by evidence, not by a promised calendar date. The detailed physical acquisition gates remain in the [Prisma specification](prisma-platform.md).

## Stage 1 — A verified current/spectroscopy simulation bridge

**First implementation delivered:** [mrs-jax PR #1](https://github.com/m9h/mrs-jax/pull/1) adds exact piecewise-constant phase integration with explicit SI units/sign, RF sensitivity and delays; phase inputs for MRSI/EPSI; native JAX multimetabolite MRSI; and per-sample EPSI timing. It also fixes the DFT grid/centering, raw-input mutation and NumPy/JAX short-basis disagreement. See the [six-control synthetic benchmark](../data/synthetic_mrs_current/README.md). The benchmark pins commit `138ffb863c1a4ffbfb1297d2b66c4121781fbff6`; it does not depend on that PR being merged. Validation: 120 mrs-jax tests passed with 26 skips in Linux CI, plus all 22 workbench tests passed locally.

The benchmark combines this repository's closed-loop field model with mrs-jax. It separates null, phase-only, concentration-only, combined, reversed-current and delayed-current conditions. Phase-aware fitting uses **known** phase and a matched illustrative basis, so this establishes software consistency rather than realistic mI quantification or blind recovery.

**Next within this stage:** sequence-specific coupled-spin mI/Cr/NAA bases; independent reference fitting; realistic linewidth, baseline/lipid, B0 drift and complex noise; joint inference with unknown current amplitude/delay. Acceptance requires null false-positive rates, concentration bias and uncertainty coverage over held-out conditions. The simple EPSI inverse FFT must be replaced or extended for sample-timed model-based reconstruction before claiming accurate dynamic-EPSI recovery.

## Stage 2 — Conductivity and complete-circuit field benchmarks

**Two tank milestones delivered:** the [3D conservative scalar FV reference](conductive-tank.md) and an independent [DOLFINx/Gmsh complete-electrode FEM implementation](fem-complete-electrode.md). The FEM model adds finite electrode patches, Ω·m² contact impedance, piecewise scalar/SPD tensor conductivity, CG1/CG2 and edge-refined CAD meshes. Analytical homogeneous/layered checks, a layered FV/FEM comparison, tensor rotation, reciprocity and energy/current checks pass. There are 48 distinct tests across the two environments. Red–green histories are retained for [FV](tdd-conductive-tank.md) and [FEM](tdd-fem.md). Complete-circuit magnetics and anatomical validation remain open.

Deliver versioned geometry, conductivity, electrodes, mesh, E/J/B fields and convergence reports. Require charge conservation, mesh convergence, polarity/orientation checks and agreement with independent solutions. Only then add conductivity inference; hold out geometries and injection patterns, and show cases where B-parallel data cannot identify the requested parameters. sbi4dwi supplies candidate microstructure priors and inference components after these forward checks.

The **next coding ticket** is conservative FEM current reconstruction with cell/interface charge checks and error estimates. Then compare the same CEM problem independently in GetDP, specify outgoing/return leads, and validate complete-circuit magnetic fields before connecting tank J to MR encoding. The [PDE/FEM roadmap](fem-roadmap.md) defines PETSc/MPI scaling, MFEM evaluation and derivative gates. Neither tank exports B yet. The FV fixed-current conductivity-scaling null remains an identifiability check for later inversion.

## Stage 3 — Measured timing and current, before the scanner

Build a deterministic controller with preloaded schedules and explicit arm/run/stop/fault states. Log scanner-marker input, delivered current and voltage on an independently observed clock; host Python handles configuration and data transfer. Exercise missed/extra triggers, delay/drift, current settling, saturation and shutdown on characterized bench loads.

Exit evidence: latency/jitter and drift distributions, waveform-integral error, calibration uncertainty and fault-injection records. Set numerical acceptance thresholds from the intended phase error budget before testing. Physical interface choice depends on the site's and Siemens' documented ports, isolation requirements and sequence timing semantics.

The [Siemens technical brief](prisma-platform.md#technical-brief-for-siemens-and-the-local-mr-team) is prepared but unsent. Needed inputs are the exact software/WIP build, sequence-development route, supported trigger/gating interfaces, raw-data timing labels and redistribution terms. No undocumented pinout or arbitrary external RF/gradient control is assumed.

## Stage 4 — Prisma phantom releases

First scan a registered closed-loop reference with current-off, reversed polarity, repeated geometry and deliberate timing offsets. Predict complex phase from the **measured** waveform and actual event times. Then progress to a conductive phantom with independently characterized conductivity and multiple injection patterns.

Exit evidence: spatial residual maps, repeatability, timing uncertainty, held-out predictions and documented failures. The local MR engineering process assesses the actual assembly/lead/coil/sequence configuration before scanning. Publish raw complex data, geometry, current/trigger logs, calibration, sequence/protocol information, reconstruction and QC under the [open-data contract](open-data.md). Synthetic releases cannot satisfy this measurement gate.

## Stage 5 — Reproduction and quantitative inference

Freeze the analysis and test a new session and independent site/configuration. Evaluate predictive residuals and uncertainty coverage rather than conductivity-map appearance. After successful phantom reproduction, plan human MRCDI and longitudinal MRSI with the site's approved equipment and study pathway. Metabolite response is a separate biological endpoint, not a direct electrical-field meter.

## Repository ownership and present status

| Repository/component | Responsibility | Present status |
|---|---|---|
| **neuroelectric-mri** | Literature, conventions, physical references, integration benchmarks, data contracts and scanner requirements | Analytical/MRSI references, conservative scalar FV and DOLFINx tensor-CEM tanks with synthetic releases and red–green tests |
| **mrs-jax** | Spectroscopy phantoms, forward models, processing and phase/timing support | Implementation on `feature/current-encoding`; [PR #1](https://github.com/m9h/mrs-jax/pull/1) open for review, not merged |
| **SimNIBS** | External reference for anatomical electrical/field models | Selected reference; no tank/head comparison executed here yet |
| **sbi4dwi** | Candidate microstructure priors, inference and future scanner/Bloch integration | Audited previously; no modifications in this work package |
| **Bench/Prisma equipment** | Independent waveform delivery/observation and scanner synchronization | Requirements drafted; no hardware built or scanner run |

This order makes every later physical measurement test an explicit prediction. It also prevents phase artifacts from being mistaken for metabolic changes while the spectroscopy and MREIT branches develop together.
