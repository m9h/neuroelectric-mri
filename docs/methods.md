# Names, mechanisms, and present evidence

Assessment date: **2026-10-23**. “SOTA” here means the strongest available methods and evidence, including negative experiments. It does not imply a deployable neural recording modality.

## Vocabulary

| Name | Measurement or proposed mechanism | Relationship to the central question |
|---|---|---|
| ncMRI, neural/neuronal-current MRI; NCI | Endogenous current → magnetic field → MR phase or dephasing | Core topic |
| msMRI, magnetic-source MRI | Historical name associated with magnitude effects of neural fields | Especially the Xiong lineage |
| Direct neuronal MRI / direct fMRI | Broad, inconsistently used description | Always identify the actual contrast mechanism |
| SIRS, stimulus-induced rotary saturation | Spin-locked magnetization responds to a resonant low-frequency field | Direct-field sensing candidate |
| REX, rotary excitation; NEMO, neuro-electro-magnetic oscillations | Related spin-lock approaches with different preparations and readouts | Do not merge all implementations or sensitivity numbers |
| ULF-NCI, dcNCI | Ultra-low/near-zero-field encoding of current-generated fields | Specialized hardware and protocols |
| LEI, Lorentz-effect imaging; NAMO | Current × magnetic field → mechanical effects encoded by gradients | Different physics from simple field-to-phase detection |
| DIANA, direct imaging of neuronal activity | Rapid repeated acquisition; originally interpreted as activity-linked relaxation changes | Original central paper retracted |
| TRIGGER | Time Resolved ImaGing of Global Electroneurophysiological Record | Related rapid human-imaging claim, preprint in this survey |
| qsMRI, quantum-sensing MRI | Time-resolved FID phase attributed to neuronal fields | Emerging 2026 preprint claim |
| MRCDI / MR-CDI | MR measurement of imposed-current fields; current-density inference | Useful calibration branch |
| MREIT / MR-EIT; DT-MREIT | Low-frequency conductivity inversion, optionally constrained by diffusion tensors | Related inverse problem, not spontaneous-current detection |
| fMREIT | Activity-dependent conductivity changes probed with an imaging current | Functional conductivity route |
| EPT / MREPT | Electrical properties inferred from RF fields at the Larmor frequency | Different frequency regime from MREIT |
| CTI, conductivity tensor imaging | Model-dependent conductivity inference combining MR measurements | Check each paper's inputs and frequency assumptions |
| dfMRI / DW-fMRI | Activity-linked diffusion/microstructure contrast | Adjacent, not a direct current measurement |

## Direct endogenous fields

Conventional phase and magnitude approaches have compelling physical foundations and positive controlled-preparation reports. The blood-free cultures and isolated cerebellum are valuable, but neither establishes routine cortical recording in humans. [Petridou et al.](https://doi.org/10.1073/pnas.0603219103), [Sundaram et al.](https://doi.org/10.1016/j.neuroimage.2016.01.059)

The forward problem must include intracellular and extracellular return currents, cancellation within a voxel, temporal synchrony, and source orientation. Realistic cell simulations show why replacing tissue by one favorably oriented current dipole can be misleading. [Cassarà et al.](https://doi.org/10.1016/j.neuroimage.2007.08.048)

Spin-lock methods are an especially concrete experimental development path. [Truong et al. (2019)](https://doi.org/10.1002/mrm.27654) reported NEMO contrast during eyes-open/closed alpha modulation in eight volunteers and a 0.06 nT phantom result. That phantom number is protocol dependent and cannot be compared directly with another study's detection threshold without matching voxel size, averaging, bandwidth and decision rule.

[Capiglioni et al. (online 2024; issue 2025)](https://doi.org/10.1002/nbm.5280) reported variable visual-stimulation SIRS results: four of twelve individuals showed significant RFR contrast, while reliability and group interaction evidence were limited. It is a pilot, not settled validation.

The [2026 MRI–MEG benchmark](https://doi.org/10.1002/hbm.70624) found no consistent in vivo activation with balanced REX/SIRS at 3 T. Phantom thresholds were approximately **0.2 nT for REX** and **0.6 nT for SIRS**, versus a **~0.07 nT MEG-derived local-field estimate**. These are study-specific limits, not fundamental bounds or direct intracortical field measurements. Its figure caption gives 0.28 nT for the last significant REX point; retain that distinction when reproducing it. [Analysis code](https://github.com/milecap/MRIMEG_pipelines)

[Gregersen et al. (2026)](https://doi.org/10.1002/mrm.70116) extend the measurable frequency range using off-resonance spin-lock preparations, tested with KomaMRI and a loop phantom. Broader frequency access trades against sensitivity. This is an acquisition advance rather than new human neuronal-current validation.

Ultra-low-field approaches can change encoding and reduce susceptibility-related complications, but have their own polarization, sensor-noise and acquisition limits. A [2020 steady-state SIRS study](https://doi.org/10.1088/1741-2552/ab87fe) demonstrated phantom detection down to 7 nT in its setup; a [2021 3D dcNCI study](https://doi.org/10.3389/fphy.2021.647376) combined simulations and phantoms. Neither result supplies a general human detection threshold.

## DIANA, TRIGGER and qsMRI

DIANA's original result should be labeled **retracted**, including in imported bibliographies. Choi et al. found no response without data preselection; Phi Van et al. produced DIANA-like signals with non-neural explanations supported by controls and simulations. [Choi et al., 2024](https://doi.org/10.1126/sciadv.adl0999), [Phi Van et al., 2024](https://doi.org/10.1126/sciadv.adl2034), [2025 retraction](https://pubmed.ncbi.nlm.nih.gov/40997196/)

TRIGGER's [2023 preprint](https://doi.org/10.1101/2023.05.31.542967) reports 18 participants, 5 mm spatial resolution and 1.4 ms temporal sampling. Images combine repeated stimulus-synchronized acquisitions. Those numbers do not mean a new complete image from one unrepeatable neural event every 1.4 ms. No peer-reviewed successor or independent confirmation was identified in this search; the DIANA retraction is not itself a retraction of TRIGGER.

[Qian et al., 2026](https://arxiv.org/abs/2601.16423v2) propose qsMRI using water-proton FID phase, with simulations, phantoms and human observations. Treat these as author claims pending independent tests. A [September 2026 mini-review](https://arxiv.org/abs/2609.21106v2) is by the same lead author, so is not replication. Our proposed checks are voxel/coil averaging, intra- versus extra-axonal cancellation, FID multi-component beating, drift, physiological motion, and independently measured electrophysiology. Calling proton spins “quantum sensors” does not by itself establish neural specificity.

## Adjacent routes

[NAMO (2006)](https://doi.org/10.1073/pnas.0605486103) reported human median-nerve effects using oscillating gradients and a Lorentz-effect interpretation. Peripheral-nerve results do not establish cerebral performance; an electromechanical model is required before folding this into a neural-field simulator.

[Nunes, Gil & Shemesh (2021)](https://doi.org/10.1016/j.neuroimage.2021.117862) provide evidence for a rapid diffusion-functional response related to neuromorphological coupling. That is highly relevant to sbi4dwi's microstructure models, but remains a separate mechanistic hypothesis from magnetic current sensing. Likewise, [functional MREPT (2024)](https://doi.org/10.1016/j.jneumeth.2024.110288) reports visual-task conductivity changes; such RF-frequency contrast cannot alone identify endogenous currents or exclude vascular mechanisms.

For practical experiments, the evidence supports using MRI with EEG/MEG or invasive recordings as independent comparators. The research opportunity is to determine which electrical parameters MRI can actually constrain, at what bandwidth and with what confounds.
