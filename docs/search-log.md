# Search and maintenance record

Initial survey: **2026-10-23**. Last targeted search: **2026-10-25**, covering active MREIT groups, software, stimulation spectroscopy, scanner synchronization and open-data requirements. This is a curated research survey, not a systematic review or citation-count analysis. “No established method” and “no independent replication identified” are bounded assessments of the material located, not proofs that no other report exists.

## 25 October software and platform audit

Added nine papers, bringing the registry to 52. The [ecosystem audit](mreit-ecosystem.md) separates published software descriptions from verified downloads/releases and distinguishes injected-current MREIT from adjacent CTI/EPT work. No private project history was added as scientific evidence.

Inspected SimNIBS 4.6.0 release metadata/tree and `calc_B` source, EPTlib metadata/README/license, ROAST metadata, DTU UQ thesis code documentation, official FSL-MRS documentation, Pulseq trigger/output source and CMRR's distribution page. These dependencies were not installed or executed. The MRCI/CoReHA download site was inaccessible. OSF GCU5B/DCWUN links were verified through papers but their file contents and licenses were not retrieved. The 2025 Copenhagen paper explicitly distinguishes public sub-study 3 code from reconstruction/optimization scripts available on request pending release.

Searches included full method names to avoid acronym collisions with real-estate investment trusts, multi-resolution EIT and unrelated MRCDI repositories. Contemporary affiliations were checked through institution records; historical collaborators are not all assumed to remain at their paper-era institution. The Prisma hardware document is a proposed requirements/verification plan; no scanner, stimulator or sequence was operated.

Follow-up: inspected m9h/mrs-jax at commit `3220889281da83f49c618cc5cc623fb6b56b3a33` after its omission from the initial software audit was identified. Retrieved source read-only and ran three upstream simulation test files from a temporary snapshot: 18 passed, 5 skipped for unavailable WAND/ISMRM data. Additional small runtime checks covered gradients, reconstruction input mutation and short-basis behavior. The [audit](mrs-jax-integration.md) distinguishes implemented prototypes from stale README roadmap entries and revises the spectroscopy development plan around mrs-jax. No mrs-jax source changes or full-suite validation were performed.

Implementation follow-up on the same date: cloned mrs-jax, followed its feature-branch/PR contribution workflow, and implemented [PR #1](https://github.com/m9h/mrs-jax/pull/1) at `138ffb863c1a4ffbfb1297d2b66c4121781fbff6`. Baseline full suite: 92 passed, 26 skipped; initial modified local suite: 119 passed, 26 skipped. Linux CI exposed a pre-existing mixed-precision MEGA-PRESS comparison; explicit single/double-precision cases now pass with precision-specific tolerances, and the final Linux suite is 120 passed, 26 skipped (JAX 0.12.0.1). An existing JAX float64-to-float32 warning remains. All 22 neuroelectric-mri workbench tests passed with the new dependency available. This is software/analytical verification; skipped reference data were not downloaded or validated, and no physical equipment was operated.

## Search families

- `neuronal current MRI`, `neural current imaging`, `ncMRI`, `NCI`, `magnetic source MRI`, `msMRI`.
- `Xiong Fox Gao San Antonio`, `Bodurka Bandettini`, `Petridou neuronal electrical activity`.
- `spin lock neuronal current`, `SIRS`, `REX`, `NEMO`, `ultra low field dcNCI`, with 2024–2026 date terms.
- `DIANA replication retraction`, `TRIGGER MRI`, `quantum sensing MRI`, including September 2026 follow-up material.
- `MREIT`, `MRCDI`, `DT-MREIT`, `functional MREIT`, `conductivity tensor imaging`, `MREPT`, paired with 2024–2026 and human-brain terms.
- `POSSUM time-varying magnetic fields`, `Drobnjak Pell Jenkinson`, and current-timing experiments.

Primary papers, publisher records, PubMed/PMC, authors' institutional repositories, and official software records were preferred. News and bibliographic aggregators were used to discover primary sources, not as sole support for the substantive conclusions. Reviews are explicitly marked as reviews in the registry.

## Important verification decisions

- DIANA's retraction was checked against PubMed's original-paper correction linkage and the 2025 retraction notice. The Science notice itself returned an access error in this session.
- qsMRI's research report and September mini-review are preprints; PubMed/PMC indexing does not convert them into peer-reviewed journal articles.
- The 2026 SMS-EPI study was checked in full at PLOS: the current flowed in an external loop around the head. An abstract-only description could easily be mistaken for injected-current brain validation.
- The 2025 spin-lock pilot appeared online in November 2024. The 2026 off-resonance paper appeared online in 2025. The registry follows issue year and records the distinction in notes.
- Some PMC full-page requests encountered a browser check. Abstracts, indexed full-text passages and author-hosted copies supplied the recorded evidence; individual registry entries state the inspection level.
- Kamei 1999 metadata and historical role were cross-checked in author-originated references; the original full text was not inspected.
- Publication dates were separated from crawl, repository deposit and PDF download dates.
- The local sbi4dwi code audit was read-only and revision-pinned. No claim is made that all its existing tests pass or that the library already supports the proposed neural simulation.

## Resources to revisit

| Resource | Access verified here | Use |
|---|---|---|
| [MRIMEG_pipelines](https://github.com/milecap/MRIMEG_pipelines) | Public repository metadata through GitHub connector | Code associated with 2026 spin-lock study; not executed here |
| [RFR_NS_pipelines](https://github.com/milecap/RFR_NS_pipelines) | Link in the pilot paper | Earlier spin-lock processing; not fetched here |
| [SMS phantom data](https://doi.org/10.17605/OSF.IO/DCWUN) | Data availability statement at PLOS | Download/access and data license still to inspect |
| [KomaMRI](https://github.com/JuliaHealth/KomaMRI.jl) | Public repository metadata through GitHub connector | Independent Bloch-simulation candidate |
| [POSSUM documentation](https://pages.fmrib.ox.ac.uk/fsldocs-4c7e2c/utilities/possum.html) | Official FSL documentation | Historical capability; local installation not tested |

## Updating the bibliography

Edit `references/sources.json`, retaining status, evidence scope, limitations and verification level. Then run:

```sh
python3 scripts/build_references.py
python3 scripts/build_references.py --check
```

The script regenerates BibTeX, CSV and the annotated reading list. It validates internal structure, not remote link availability or scientific truth. Link-checking and publication-status updates remain deliberate research tasks. Add a retraction note wherever a retracted paper is exported; never cite only the original title as positive evidence.

For the next literature pass, prioritize independent qsMRI/TRIGGER tests, human spin-lock replication, fMREIT translation, held-out-montage conductivity validation, and release of raw complex data. A true bibliometric “peak” claim would require a separate reproducible query with acronym disambiguation and year-by-year counts.
