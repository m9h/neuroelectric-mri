# Neuroelectric MRI

**A research map and simulation workbench for neuronal-current MRI, MR current-density imaging, and electrical conductivity imaging.** Literature checked through **23 October 2026**.

The central assessment is that endogenous neuronal-current MRI remains an experimental research problem. The reviewed literature does not establish a generally reproducible human method with both useful spatial localization and millisecond neural specificity. Measuring externally imposed current fields is substantially further along. These conclusions concern different measurement problems, so this repository keeps their evidence separate.

The San Antonio paper you probably remember is **Jinhu Xiong, Peter T. Fox, and Jia-Hong Gao (2003), “Directly mapping magnetic field effects of neuronal activity by magnetic resonance imaging.”** All three were at the University of Texas Health Science Center at San Antonio. It introduced the influential **magnetic-source MRI / msMRI** claim. Xiong is the likely name; Gao is another plausible recollection. It was an early prominent report, rather than the first proposal or experiment. [Paper](https://doi.org/10.1002/hbm.10124)

Several developments change the assessment today:

- **DIANA:** the original 2022 Science paper was **retracted on 25 September 2025**. [Retraction](https://pubmed.ncbi.nlm.nih.gov/40997196/)
- **Spin-lock MRI:** a 2026 MRI–MEG study provides a useful negative human benchmark and quantitative sensitivity estimates. [Capiglioni et al.](https://doi.org/10.1002/hbm.70624)
- **qsMRI:** a 2026 preprint makes new human detection claims; independent validation remains an open question in this survey. [Qian et al.](https://arxiv.org/abs/2601.16423v2)
- **MREIT/MRCDI:** human field measurements, conductivity-model fitting, and acquisition improvements continued well beyond the 2000s. [2024 model calibration](https://doi.org/10.1162/imag_a_00176), [2026 SMS acquisition study](https://doi.org/10.1371/journal.pone.0341731)
- **POSSUM is directly relevant:** its time-varying-field extension already explored current-induced MRI effects. [Drobnjak, Pell & Jenkinson, 2010](https://doi.org/10.1016/j.mri.2010.03.029)

## Read the repository

| File | Purpose |
|---|---|
| [History](docs/history.md) | San Antonio, Bandettini, early claims and later developments |
| [Methods and present evidence](docs/methods.md) | Names, mechanisms, evidence boundaries, current approaches |
| [MREIT and conductivity](docs/mreit.md) | What survived, functional MREIT, DT-MREIT and EPT distinctions |
| [Simulation roadmap](docs/simulation-roadmap.md) | Inspected sbi4dwi code, missing physics, concrete benchmarks |
| [Annotated reading list](docs/reading-list.md) | Study-level findings and limitations |
| [Search record](docs/search-log.md) | Search scope, dates, access limits and update procedure |
| [Source registry](references/sources.json) | Machine-readable bibliography and evidence annotations |
| [BibTeX](references/references.bib) / [CSV](references/literature.csv) | Citation-manager and spreadsheet exports |
| [Phase-budget example](examples/phase_budget.py) | Executable sanity checks for field-to-phase sensitivity |

## Why this name?

`neuroelectric-mri` is a useful umbrella: it accommodates endogenous current detection, injected-current calibration, and tissue conductivity without declaring any contested technique successful. The conventional narrow literature term is **neuronal-current MRI (ncMRI)**. Alternative repository names are `neuronal-current-mri` for a narrower survey or `bioelectric-mri` for eventual cardiac and peripheral-nerve coverage. NEMO, NAMO, and DIANA already identify particular methods and would make confusing umbrella names.

## Run the small computational example

Python 3.10+; standard library only:

```sh
python3 examples/phase_budget.py
python3 examples/phase_budget.py --csv data/phase_budget.csv
python3 -m unittest discover -s tests -v
python3 scripts/build_references.py --check
```

The example computes ideal phase budgets and exact cancellation cases. It is **not** a neuronal-current simulator, an MREIT reconstruction, or a reproduction of a published experiment. The larger simulator is a proposed roadmap. No changes to sbi4dwi are required to run this repository.

For a research program, start with calibrated current phantoms and artifact controls; then test whether plausible neural sources remain identifiable after voxel averaging, acquisition, and nuisance variability. Simulation-based inference becomes useful once this forward problem and its null cases are credible.
