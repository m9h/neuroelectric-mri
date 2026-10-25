# Open MREIT data contract — draft 0.1

Scope: a proposed phantom-first release format, drafted 25 October 2026. No acquired dataset has been released by this repository. The checked-in analytical examples are **synthetic reference calculations**, not scanner data or a conductivity ground truth.

## Release contents

| Asset | Required information |
|---|---|
| Raw MR data | Vendor raw data where redistribution is permitted; complex multicoil samples, sampling trajectory, noise/reference scans, reconstruction prerequisites and acquisition flags; MRD conversion when supported |
| Sequence/protocol | Sequence or protocol file, checksum, software/interpreter/WIP versions, RF/gradient/ADC timing, all echoes, current-encoding windows, slice/shot ordering, dummy acquisitions, hardware raster times |
| Current waveform | Separate commanded and independently measured current; voltage/compliance state; SI units; sampling rate, bandwidth/filter delay, calibration and uncertainty |
| Synchronization | Raw trigger edge timestamps, event identity/polarity, clock origin and mapping, drift/jitter estimates, missing/extra edges, aborts and run state |
| Geometry | Phantom CAD, electrode and lead/return paths, material assignments, fiducials, image/mesh/scanner transforms, axis handedness, B0 direction, geometry uncertainty |
| Conductivity reference | Independently measured reference values with frequency, temperature, method, calibration, uncertainty and repeatability; avoid calling nominal saline concentration ground truth |
| Derived products | Complex images, phase differences with sign convention, B-parallel maps and uncertainty, masks, E/J predictions, reconstructed conductivity, residuals and QC |
| Provenance | Acquisition dates, hardware/firmware and analysis revisions, dependency lockfiles/container digest, parameter files, seeds for synthetic data, file hashes and dataset version |
| Controls | Current off, reverse polarity, repeat runs, timing offsets, orientation/geometry nulls, instrument faults and excluded runs with reasons |
| Reuse | Dataset license, citation/DOI, data dictionary, acquisition instructions, exact reproduction command and limitations; code license handled separately |

The raw record must distinguish **commanded**, **measured**, **simulated** and **reconstructed** quantities. A current-source setting is not a current measurement; a SimNIBS prediction is not experimental ground truth.

## Coordinate and time contracts

Use explicit SI units internally. Store the transform from geometry coordinates to scanner coordinates and to NIfTI world coordinates. Define B0's unit vector in the same frame as the predicted field before projecting it. Record the Fourier/complex-signal sign convention and whether a difference means `phase(+I) - phase(-I)` or the opposite.

Keep the device's raw monotonic timestamps. Store a calibrated mapping to the scanner-event timebase with uncertainty, including drift rather than only a single offset. Event zero must be defined explicitly; retained-image zero and first hardware trigger can differ because of dummy scans. For each encoded measurement, preserve the waveform and signed sequence-sensitivity interval that generated it. Phase bias from measured timing differences is a first-class QC output.

## Directory proposal

```text
dataset/
  README.md
  dataset_description.json
  manifest.json                 # schema/version, provenance, checksums
  sourcedata/
    scanner_raw/
    sequence/
    hardware_logs/
    geometry/
    conductivity_calibration/
  derivatives/
    reconstruction/
    current_field/
    conductivity/
    qc/
  code/
```

This layout is an initial release convention, **not a claim of BIDS validation**. Use valid BIDS MRI representations for eligible converted data and documented derivatives. Preserve additional instrument records in sourcedata with explicit sidecars rather than inventing standardized suffixes. Track [BEP037 for non-invasive brain stimulation](https://bids.neuroimaging.io/extensions/beps/bep_037.html); the official page describes its March 2026 version as work in progress. Use [MRD/ISMRMRD](https://ismrmrd.github.io/) for portable raw MR representation when converters preserve the necessary information. Keep original raw data and test conversions for loss of phase, labels and timestamps.

## Release acceptance

1. A fresh environment reproduces at least one documented field/phase comparison from the release inputs.
2. File hashes, units, geometry transforms, timebase mappings and polarity labels validate; missing records are explicit.
3. QC includes a conservation/forward-model check where applicable, timing failures, current-monitor saturation and phase residuals. Failed runs remain discoverable.
4. Train/test partitions separate physical phantoms, sessions or sites; slices from the same acquisition are not treated as independent validation subjects.
5. Synthetic and measured releases have distinct identifiers. A simulated benchmark cannot satisfy a measured-data milestone.
6. Redistribution rights for raw data and scanner-related assets are recorded. Any future human release uses consent and privacy handling designed for sharing, including raw headers and potentially identifying anatomy.

Host source and small fixtures on GitHub; archive substantial immutable releases on a research data service with a DOI, checksums and a versioned dataset card. OSF GCU5B/DCWUN are starting points for external comparisons (see [the audit](mreit-ecosystem.md)), not data already mirrored in this repository.
