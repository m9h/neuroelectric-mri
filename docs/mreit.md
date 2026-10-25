# MREIT: a continuing field with a more specific role

See the [25 October 2026 groups/software audit](mreit-ecosystem.md), [Prisma platform plan](prisma-platform.md) and [open-data contract](open-data.md) for the practical next steps.

The 2005–2012 period was a major development era, but this survey does **not** establish a bibliometric peak. MREIT did not disappear. Its current relevance includes stimulation dosimetry, low-frequency tissue properties, and validation of volume-conductor models. Human current-field measurements have stronger practical evidence than spontaneous neuronal-current imaging.

## What the scanner measures

An imposed current pattern produces a magnetic field. With a fixed head orientation, conventional MRI is primarily sensitive to the perturbation parallel to the main field, denoted Bz. For a constant field over an effective encoding time, opposite current polarities give a phase difference of magnitude

```text
Δφ = 2 γ Bz T_eff
γ = 2π × 42.57747892e6 rad s⁻¹ T⁻¹
```

The sign depends on phase and polarity conventions. General waveforms require integration with the sequence's signed sensitivity function; T_eff is not automatically TE. Reversal suppresses common background phase, but does not automatically cancel drift between acquisitions, motion, lead fields, or stimulus-linked interference. The double-echo EPI implementation in [Gregersen et al. (2024)](https://doi.org/10.1162/imag_a_00176) provides a concrete acquisition/model-fitting example.

For a quasistatic, passive conductor with conductivity tensor σ,

```text
J = −σ ∇u
∇·J = 0
∇×B = μ₀ J
```

Boundary conditions specify electrode currents and insulation. MREIT seeks σ; MRCDI seeks information about J. The observed Bz is not itself the full current-density vector. Recovering J from a curl requires unavailable components or extra assumptions; conductivity inversion adds further boundary, geometry and parameter-identifiability issues. For active tissue, an impressed source or a bidomain formulation must be added.

## Why translation was difficult

The original developers' [2008 review](https://doi.org/10.1088/0967-3334/29/10/R01) and [2011 mathematical account](https://doi.org/10.1137/080742932) explain the progression toward Bz-based methods that avoid rotating the subject. The continuing obstacles are small phase effects at practical currents, spatial derivatives that amplify noise, incomplete field observations, electrode/lead modeling, and long acquisitions. Detailed anatomical images do not remove those inverse-problem limitations.

[Kim et al. (2009)](https://doi.org/10.1109/TMI.2009.2018112) reported human-leg conductivity images with 1.7 mm pixels in four volunteers. Their protocol used brief currents up to 9 mA in the leg; this is historical experimental context, not a transferable brain-stimulation setting. [Oh et al. (2004)](https://doi.org/10.1002/mrm.20091) is a useful earlier 3 T milestone.

## What has progressed since then?

| Study | Contribution | Boundary of the evidence |
|---|---|---|
| [Göksu et al., 2018](https://doi.org/10.1016/j.neuroimage.2017.12.075) | Human in vivo brain MRCDI | Measures imposed-current fields; current reconstruction is model dependent. |
| [Chauhan et al., 2018](https://doi.org/10.1109/TMI.2017.2783348) | Human-head DT-MREIT | Tensor inference uses diffusion-based assumptions. |
| [Göksu et al., 2021](https://doi.org/10.1002/mrm.28944) | Reported 84 pT field sensitivity in 4:20 min | An acquisition-specific sensitivity, not millisecond neural recording. |
| [Göksu et al., 2023](https://doi.org/10.1002/mrm.29780) | Volumetric/SMS acquisitions; 67 pT reported for three distant slices, 2 × 2 × 3 mm³, 10 min | Sparse coverage and averaging matter. |
| [Gregersen et al., 2024](https://doi.org/10.1162/imag_a_00176) | Double-echo EPI and individualized conductivity optimization | Five-subject study; fitting does not uniquely identify every tissue parameter. |
| [Mosayebi-Samani et al., 2025](https://doi.org/10.1162/imag_a_00481) | Tests anisotropy and reconstruction assumptions | Tissue conductivity uncertainty mattered more than DTI anisotropy in the tested cortical tES settings. |
| [Cunha et al., 2026](https://doi.org/10.1371/journal.pone.0341731) | SMS-EPI with MB factors up to 5 and adequate slice gaps | Phantom and five-human-head **external cable-loop** experiments; not endogenous currents, and not a new through-head tES validation. |

The 2026 SMS study publishes phantom data at [OSF](https://doi.org/10.17605/OSF.IO/DCWUN); human data are restricted. Its ≥12 mm interslice gap condition matters: “more coverage” is not contiguous whole-brain coverage. The [2026 conductivity-imaging review](https://doi.org/10.1002/jmri.70279) places this work alongside EPT and clinical research applications, while identifying outstanding translation problems.

## Functional MREIT is an important bridge

fMREIT probes activity-linked **conductivity changes** with an applied current. Changes in membrane conductance can alter the path of that probing current and hence the MR phase. This differs from passively observing the magnetic field of native neural currents.

[Sadleir et al. (2017)](https://doi.org/10.1016/j.neuroimage.2017.08.004) reported activity-related MREIT phase variability in isolated Aplysia ganglia, including current-off and solution controls. It is an in vitro result, not evidence for whole-brain human functional MREIT. Their [2019 active-bidomain simulation](https://doi.org/10.1002/mrm.27351) attributes the modeled contrast primarily to membrane-conductivity changes rather than source strength, and notes that the probing current can itself affect neural activity. This is a particularly relevant paper for a mechanistic simulator.

## The sbi4dwi connection: D is not σ

Water diffusion and electrical conductivity can share microstructural constraints without being interchangeable quantities. Diffusion MRI measures water displacement; conductivity depends on mobile charge, membranes, frequency and compartments. DT-MREIT commonly constrains tensor geometry using DTI and estimates scale factors from electrical measurements. It does not directly measure an unconstrained six-component conductivity tensor in each voxel.

[Hakhu et al. (2025)](https://doi.org/10.3389/fradi.2025.1492479) compare conductivity estimates derived from NODDI, SMT and SANDI and find substantial model dependence. Some comparisons also involve different datasets; they are not a pure same-subject model bake-off. The [2025 anisotropy study](https://doi.org/10.1162/imag_a_00481) further cautions that realistic Bz-derived reconstructions are less favorable than reconstructions supplied with ideal current densities.

Our proposed inference target is therefore a **joint posterior over microstructure, conductivity mapping, electrode/contact parameters and measurement nuisance terms**, with held-out current patterns. Do not report conductivity uncertainty as though choosing a diffusion model had removed it.

## MREIT versus EPT and CTI

MREIT probes low-frequency electrical behavior using external currents. EPT estimates RF-frequency properties from transmit/receive field information: approximately 128 MHz at 3 T and 298 MHz at 7 T for protons. CTI methods may combine diffusion and RF conductivity with compartment models to infer lower-frequency properties. Frequency conversion is a modeling assumption, not an identity. [Conductivity-model comparison](https://doi.org/10.3389/fradi.2025.1492479)

Functional EPT reports, such as [Kim et al. (2024)](https://doi.org/10.1016/j.jneumeth.2024.110288), belong in the watch list, but do not demonstrate that millisecond membrane currents have been directly imaged. Conventional electrode EIT, MR-EPT, and magnetic resonance elastography (MRE) are also distinct modalities despite overlapping abbreviations.
