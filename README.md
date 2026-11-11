# Dementia prevalence in Brazil: why national estimates disagree

Reproducibility repository for:

> Pereira da Silva AM, Santos DH, Lima Júnior DVS, Tudella GCN, Cardoso LJ, Ferreira GD, de Castro WCN, Pinto ELO, Lee AM, Ahmed A, Perry G, Duarte FS. **Dementia prevalence in Brazil and why national estimates disagree: a systematic review and meta-analysis with bias analysis.** Manuscript under review.

PROSPERO registration: CRD420261341673.

The repository contains the curated study-level data, every analysis script, and the outputs behind each table and figure of the article and its Supplementary material. A single command reruns the whole analysis.

## Contents

- [Quick start](#quick-start)
- [Data](#data)
- [Repository structure](#repository-structure)
- [Analysis pipeline](#analysis-pipeline)
- [From outputs to the article](#from-outputs-to-the-article)
- [Computational environment](#computational-environment)
- [Reproducibility checks](#reproducibility-checks)
- [Licence and citation](#licence-and-citation)

## Quick start

```bash
git clone https://github.com/Andersonmatheuss/dementia-prevalence-brazil.git
cd dementia-prevalence-brazil
Rscript -e 'install.packages("renv"); renv::restore(lockfile = "code/renv.lock")'
Rscript -e 'install.packages(c("sf", "geobr", "patchwork", "jsonlite"))'
python3 -m pip install -r code/requirements.txt

code/run_all.sh              # one pass: output/tables, output/figures, supplementary/supplement_tables.md
code/run_all.sh --check      # two passes; every output must be byte-identical between them
```

`PY=/path/to/python code/run_all.sh` selects the Python interpreter. The map figure downloads Brazilian boundaries through `geobr` and needs internet access.

## Data

All inputs are **published aggregate data**:

- study-level estimates and age-, sex-, schooling- and class-specific counts extracted from the full texts of the included reports, each value with the table and page it came from (`data/raw/extraction/fulltext_strata.csv`, `fulltext_study_estimates.csv`);
- ELSI-Brazil, entered as the published weighted estimates (Bertola et al., 2023), never as microdata;
- population counts from the 2022 Census (IBGE SIDRA) and the 2024 revision of the IBGE population projections.

No individual-level data are included. Full texts of the included reports are not redistributed for copyright reasons; each report is identified by PMID or DOI in `data/processed/extraction.csv`. The search strategies are given in eMethods 1 of the Supplementary material.

| File | Content |
|---|---|
| `data/processed/extraction.csv` | One row per report: population, design, ascertainment, published estimates, role in the analysis |
| `data/processed/analysis_sets.csv` | Membership of each independent population in the primary and sensitivity sets |
| `data/processed/risk_of_bias.csv` | Risk-of-bias ratings used by an earlier version of the analysis |
| `data/raw/extraction/` | Full-text extraction, adjudication of discrepancies, JBI risk of bias, checks against the sources |
| `data/raw/ibge_projection_2024/` | IBGE population projections, 2024 revision |
| `data/raw/Dementia_prevalence/ibge_sidra_inputs/` | 2022 Census and earlier projection tables by age group |
| `data/raw/Dementia_prevalence/dementia_outputs/`, `data/raw/Prevalence_Dementia_23_03/`, `data/raw/dementia-prevalence-brazil-addc1b1/` | Inputs and code of an earlier version of the analysis, kept to verify the data chain (scripts `02a`–`02c`) |
| `data/data_dictionary.csv` | Variables of the processed files |

## Repository structure

```
code/
  run_all.sh            single entry point; runs every step in order
  helpers/setup.R       seed, paths, analysis sets and shared functions
  helpers/figure_style.py
  renv.lock             R 4.6.1 and exact package versions
  requirements.txt      Python dependencies
data/
  processed/            curated analysis data
  raw/                  extraction and population inputs (read-only)
output/
  intermediate/         derived inputs and cached Bayesian fits (models/)
  tables/               results tables, one CSV per analysis
  figures/              figures produced by the pipeline
supplementary/          created by the pipeline (supplement_tables.md)
```

## Analysis pipeline

`code/run_all.sh` runs the steps below in this order; each step reads only `data/` and files written by earlier steps.

| Step | Scripts | What it does |
|---|---|---|
| Data chain | `02a`–`02c` | Checks the raw data chain and reproduces an earlier version of the analysis (frequentist and Bayesian) |
| Core analyses | `03a`, `03c`, `03f`, `03g`, `03j` | Multilevel logit model of age-specific prevalence; analysis sets; subgroups; authors' estimates as primary input |
| Population | `03b`, `03m` | 2024 IBGE projections; burden by macroregion and federative unit |
| Case ascertainment | `03d`, `03e`, `03s` | External comparison with ELSI-Brazil; deterministic NPV scenarios and probabilistic bias analysis; PPV of unexamined screen-positive participants |
| Bayesian models | `03i`, `03l`, `03r` | Bayesian age model, leave-one-population-out and prior sensitivity |
| Further analyses | `03h`, `03k`, `03n`, `03o`, `03p`, `03q` | Value of information; Pietà out-of-sample check; trend scenarios; reconciliation of national estimates; literacy standardisation; age-band contributions |
| Tables and figures | `04a`, `05a`, `05c`, `05d` | Supplementary tables (`supplementary/supplement_tables.md`) and figures |
| Checks | `06a`, `06b` | Reproduction checks; underdiagnosis by national denominator |

The project seed is `20261019`. Stochastic steps (bias analysis, simulations, Bayesian fits) set their own seeds before sampling.

The Bayesian fits (`02c`, `03i`, `03l`) reuse the models cached in `output/intermediate/models/` while data and model code are unchanged (`brms.file_refit = "on_change"`). Delete that folder to refit from scratch; this takes considerably longer.

## From outputs to the article

| Article item | Source in `output/` |
|---|---|
| Table 1 | `data/processed/extraction.csv`, `data/processed/analysis_sets.csv` |
| Table 2 | `tables/primary_age_strata.csv`, `author_best_pooled.csv`, `author_best_burden_projections.csv`, `author_best_projections_regions.csv`, `pba_screening_correction.csv`, `primary_burden_estimands.csv` |
| Figure 2 | `tables/primary_age_strata.csv`, `corrected_subgroups.csv`, `corrected_subgroup_tests.csv` |
| Figure 3 | sources of Table S13; panel C, `tables/primary_voi.csv` |
| Figure 4 | `tables/reconciliation_bridge_2022.csv` |
| Table S8 | `intermediate/author_best_estimates_s6.csv`, `data/raw/extraction/fulltext_study_estimates.csv` |
| Table S9 | `tables/author_best_pooled.csv`, `author_best_sensitivity_sets.csv`, `corrected_overall_with_epifloripa.csv` |
| Table S10, Figure S7 | `tables/corrected_overall_by_set.csv`, `corrected_metaregression_diagnostic.csv`, `corrected_age_strata*.csv`, `corrected_national_burden_2022.csv` |
| Table S11 | `tables/corrected_subgroups.csv`, `corrected_subgroup_tests.csv` |
| Table S12, Figure S6 | `tables/primary_bayesian.csv`, `bayesian_corrected_models.csv`, `bayesian_leave_one_study_out.csv`, `primary_loo_prior_sensitivity.csv`, `intermediate/diagnostics/` |
| Table S13 | `tables/primary_burden_estimands.csv`, `primary_burden_sensitivity_sets.csv`, `author_best_burden_projections.csv`, `triangulation_national_estimates.csv` |
| Table S14 | `tables/reconciliation_bridge_2022.csv` |
| Table S15 | `tables/qba_two_phase_documented_study_level.csv`, `qba_two_phase_with_epifloripa.csv` |
| Table S16 | `tables/pba_screening_correction.csv`, `pba_study_corrections.csv`, `pba_ppv_unexamined.csv`, `uncertainty_decomposition_2022.csv` |
| Table S17 | `tables/qba_elsi_algorithm_misclassification.csv` |
| Table S18, Figure S3 | `tables/primary_external_comparison_elsi.csv`, `primary_elsi_algorithms_65plus.csv`, `pieta_out_of_sample.csv` |
| Table S19 | `tables/author_best_projections_regions.csv` |
| Table S20, Figure S4 | `tables/burden_states.csv` |
| Table S21 | `tables/author_best_burden_projections.csv`, `projection_trend_scenarios.csv` |
| Table S22 | `tables/primary_voi.csv` |
| Table S23 | `tables/underdiagnosis_by_denominator_2022.csv` |
| Table S24 | `tables/education_standardised_burden_2022.csv` |
| Table S25 | `tables/age_scale_correction_strata.csv`, `age_scale_correction_burden.csv` |
| Figure S2 | `tables/pba_study_corrections.csv`, `qba_elsi_algorithm_misclassification.csv`, `underdiagnosis_by_denominator_2022.csv` |
| Figure S5 | `tables/author_best_projections_regions.csv`, `burden_states.csv`, `projection_trend_scenarios.csv` |

Tables S4 to S7 come from `data/processed/` and `data/raw/extraction/`. The flow diagram (Figure 1A) and the final layout of the article figures were prepared by the authors from these outputs.

## Computational environment

- R 4.6.1; package versions are pinned in `code/renv.lock` (among them metafor 5.0-1, brms 2.23.0, rstan 2.32.7, dplyr 1.2.1, ggplot2 4.0.3). The figure script `05d` also uses `sf`, `geobr` and `patchwork`, and `jsonlite` is used for metadata.
- Python 3 with matplotlib 3.11.1 (`code/requirements.txt`).
- Comments inside the scripts are partly in Portuguese; file names, variable names and outputs are in English.

## Reproducibility checks

- `code/run_all.sh --check` runs the chain twice and fails if any output differs between the two passes, or if code or inputs change during the run.
- `06a_check_reproduction.R` runs consistency checks on the output tables without refitting the models.

## Licence and citation

Code is released under the MIT licence (`LICENSE`). Data files in `data/` and `output/` are released under CC BY 4.0. Please cite the article (see `CITATION.cff`).

Contact: Anderson Matheus Pereira da Silva, Federal University of Pernambuco, anderson.matheuss@ufpe.br (ORCID 0009-0003-4789-6019).
