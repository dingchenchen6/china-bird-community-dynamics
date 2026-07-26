# Reproducibility guide

This document explains what is needed to re-run the analysis end to end, in what order, and what to expect at each stage. It is written so that someone who has never seen the project can follow it.

---

## 1. What is and is not in this repository

| Included | Excluded (and why) |
|---|---|
| All R pipeline code (`code/`) | Raw occurrence records — subject to provider terms (CBRC, GBIF, eBird) |
| Derived result tables (`results/`, CSV) | Posterior sample arrays and model fit objects (`*.rds`, tens of GB) |
| Manuscripts and reports | Environmental raster inputs (WorldClim, CLCD, HFI; obtain from source) |
| Main-line figures (`figures/`) | Log files and scratch output |

Because the large intermediate objects are excluded, **the postprocessing stages cannot run until the models are refitted** (stage 04). The result tables in `results/` let you audit every number in the manuscripts without rerunning anything.

---

## 2. Software environment

- **R** ≥ 4.2
- Core packages: `spOccupancy` (occupancy models), `brms` + `cmdstanr` (driver regressions), `sf`, `terra`, `exactextractr` (spatial), `vegan`, `ape`, `picante` (diversity & phylogeny), `randomForest`, `tidyverse`, `ggplot2`, `scico`, `patchwork`.
- `cmdstanr` requires a working CmdStan installation (`cmdstanr::install_cmdstan()`).
- Analyses were run on Linux (compute server) for stages 04–05 and macOS for figures and manuscripts.

Install the main dependencies:

```r
install.packages(c("spOccupancy", "brms", "sf", "terra", "exactextractr",
                   "vegan", "ape", "picante", "randomForest", "tidyverse",
                   "scico", "patchwork", "MASS", "Kendall"))
# cmdstanr from its own repo
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

---

## 3. Input data you must obtain yourself

| Input | Source | Used by |
|---|---|---|
| China Birdwatching Records Center records | https://www.birdreport.cn (terms apply) | `01_merge_birdwatch_ebird.R` |
| eBird / GBIF China records | https://ebird.org/data/download, https://www.gbif.org | `01_merge_birdwatch_ebird.R` |
| WorldClim 2.1 bioclimatic variables | https://worldclim.org | `03_prepare_environment.R` |
| CLCD China land cover (30 m) | Yang & Huang 2021, ESSD | `03d_prepare_landuse_change.R` |
| Human Footprint (annual, 2000–2018) | Mu et al. 2022, Sci Data | `03e_prepare_hfi_change.R` |
| AVONET traits | Tobias et al. 2022, Ecol Lett | `03b_extend_traits.R` |
| EltonTraits 1.0 | Wilman et al. 2014, Ecology | `03b_extend_traits.R` |
| Bird phylogeny | Jetz et al. 2012 / BirdTree | posterior propagation stage |

Place them where `code/code_v3/00_config.R` expects, or edit the path constants there.

---

## 4. Configuration

All paths, run labels and MCMC settings are centralised:

- `code/code_v3/00_config.R` — grid size, run labels (`RUN_LABEL`), MCMC settings (`FULL_N_BURN`, `FULL_N_THIN`, `FULL_N_CHAINS = 4`), seeds.
- `code/code_v3/utils_paths.R` — path resolution (`v3_file()`, `ensure_v3_dirs()`).

Set the code directory before running any script:

```bash
export V3_CODE_DIR="/path/to/china-bird-community-dynamics/code/code_v3"
export V3_RUN_LABEL="v3_full_200sp_ar1_spatial"
export V3_OMP_THREADS=8
```

---

## 5. Pipeline order

Run in this order. Stages 04–05 are the expensive ones.

| # | Script | Purpose | Typical cost |
|---|---|---|---|
| 01 | `01_merge_birdwatch_ebird.R` | Merge CBRC + GBIF/eBird, source-aware deduplication | ~1 h |
| 02 | `02_build_survey_history.R` | Build detection/non-detection survey histories on the 100-km grid × 5 periods | ~1 h |
| 03 | `03_prepare_environment.R`, `03b`–`03e` | Environmental, trait, climate-change, land-use and human-pressure covariates | ~2 h |
| 04 | `04_run_stMsPGOcc_main.R` | **Primary model** — 200 species, spatial NNGP + AR(1), 4 chains | days, high RAM |
| 04 | `04_run_tMsPGOcc_500sp.R` | Breadth extension — 500 species, temporal (non-spatial) | days, high RAM |
| 04c | `04c_combine_chains.R` | Combine per-chain fits, compute R̂ / ESS | ~1 h |
| 05 | `05_postprocess_diversity_extended.R` | Propagate occupancy posterior into all diversity metrics | hours |
| 05b/c | `05b_mcmc_diagnostic_plots.R`, `05c_ppc_bayesian_pvalue.R` | Convergence diagnostics, posterior predictive checks | ~1 h |
| 05e/f | `05e_export_trend_draws.R`, `05f_functional_trend_pdecline.R` | Theil–Sen trend draws; functional-trend posterior P(decline) | ~1 h |
| 09/14 | `09_extended_analyses_extended.R`, `14_species_trait_regression.R` | Driver variance partitioning, RF importance, trait regressions | ~2 h |
| 12 | `12_homogenization_spatiotemporal_maps.R` | Spatial beta diversity / homogenization trend | ~1 h |
| 15–17 | `15_sensitivity_3yr_window.R`, `15b`, `16`, `17` | Sensitivity: time window, breeding season, grid size, threshold | hours |
| 06/21–24 | `06_figures_publication.R`, `21`–`24_*` | Publication figures | ~1 h |

**Hardware note.** The 200-species spatial model and the 500-species temporal model were run on a compute server with large RAM; posterior arrays reach tens of GB. They are not expected to run on a laptop.

---

## 6. Random seeds and determinism

Seeds are set centrally in `00_config.R` and in `utils_seeds.R`. MCMC uses 4 chains with fixed seeds; posterior summaries are therefore reproducible up to sampler-level parallel scheduling. Re-running stage 04 will not reproduce draws bit-for-bit across different `spOccupancy` or BLAS versions, but posterior summaries should agree within Monte Carlo error.

---

## 7. Verifying results without rerunning

Every number quoted in the manuscripts can be traced to a table in `results/results_v3/`. See `results/README.md` for the data dictionary. For example:

- corrected richness per period → `table_diversity_summary_v3_full_200sp_ar1_spatial_extended.csv`
- convergence → `table_convergence_diagnostics_v3_full_200sp_ar1_spatial_4chain.csv`
- turnover/nestedness → `table_baselga_global_v3_full_200sp_ar1_spatial_extended.csv`
- effort-confound controls → `table_effort_confound_controls_v3_full_200sp_ar1_spatial.csv`

---

## 8. Known open items

Tracked in `reports/RESEARCH_SYSTEMATIC_AUDIT_20260602.md` and `reports/SERVER_RERUN_GUIDE_20260602.md`:

1. Spatial beta-diversity (homogenization) trend requires rerunning stage 12 with 4-D posterior arrays.
2. Posterior P(decline) for functional metrics — run `05f_functional_trend_pdecline.R`.
3. Data-source detection term — run `04d_source_detection_refit.R` (needs two settings aligned with stage 04).
