# Code

R pipeline for the analysis. `code_v3/` is the **current main pipeline**; `code_v4/` is an audit/repair branch; `code_v2/` is retained for version comparison (the project deliberately keeps old versions side by side rather than overwriting them).

## Conventions

- Configuration is centralised in `code_v3/00_config.R` (grid size, run labels, MCMC settings, seeds).
- Paths resolve through `utils_paths.R` (`v3_file()`, `ensure_v3_dirs()`); set `V3_CODE_DIR` before running.
- Scripts are numbered in execution order. Scripts sharing a number prefix (e.g. `03b`–`03e`) are variants/extensions of that stage.
- Most scripts carry a bilingual header stating the scientific question, inputs, workflow and expected outputs.

## Stage map (`code_v3/`)

| Stage | Script | Purpose |
|---|---|---|
| 00 | `00_config.R` | Global configuration, run labels, MCMC settings, seeds |
| 01 | `01_merge_birdwatch_ebird.R` | Merge CBRC + GBIF/eBird records, source-aware deduplication |
| 02 | `02_build_survey_history.R` | Detection/non-detection survey histories (grid × period × repeat) |
| 03 | `03_prepare_environment.R` | Grid environmental covariates (climate, topography, land cover, HFI) |
| 03b–03e | `03b_extend_traits.R`, `03c_prepare_climate_change.R`, `03d_prepare_landuse_change.R`, `03e_prepare_hfi_change.R` | Traits and change-covariates |
| 04 | `04_run_stMsPGOcc_main.R` | **Primary model**: 200 species, spatial NNGP + AR(1), 4 chains |
| 04 | `04_run_tMsPGOcc_500sp.R` | Breadth extension: 500 species, temporal (non-spatial) |
| 04b–04d | `04b_*`, `04c_combine_chains.R`, `04d_source_detection_refit.R` | Thinning, chain combination + R̂/ESS, source-detection refit |
| 05 | `05_postprocess_diversity_extended.R` | Posterior propagation into all diversity metrics |
| 05b–05f | `05b_mcmc_diagnostic_plots.R`, `05c_ppc_bayesian_pvalue.R`, `05d_compute_stability.R`, `05e_export_trend_draws.R`, `05f_functional_trend_pdecline.R` | Diagnostics, PPC, stability, trend draws, functional-trend P(decline) |
| 06 | `06_figures_publication.R` (+ `06b`–`06d`) | Core figures |
| 09 | `09_extended_analyses_extended.R` | Variance partitioning, random-forest importance, brms driver models |
| 12 | `12_homogenization_spatiotemporal_maps.R` | Spatial beta diversity / homogenization trend |
| 14 | `14_species_trait_regression.R` | Species trait–trend regressions |
| 15–17 | `15_sensitivity_3yr_window.R`, `15b_sensitivity_breeding_season.R`, `16_sensitivity_grid_size.R`, `17_sensitivity_eps_threshold.R` | Sensitivity analyses |
| 19–24 | `19_winner_loser_species_v3.R`, `20`–`24_*` | Winner/loser analysis, concept diagrams, publication figure sets |

## Shared utilities

`utils_core.R` (I/O, logging), `utils_paths.R` (paths), `utils_diversity.R` / `utils_diversity_extended.R` (diversity metrics, `theil_sen_slope`), `utils_diagnostics.R` (R̂/ESS), `utils_spatial.R` (grids, distances), `utils_mapping.R` (map theme and layers), `utils_plots_advanced.R` (raincloud/beeswarm plots), `utils_importance.R` (random-forest importance), `utils_seeds.R` (seeds).

## Mapping standards

All maps use `utils_mapping.R` and follow fixed rules: unified bbox `coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE)`, no inset/locator maps, no `NA` facet panels, missing grids shown explicitly as grey rather than dropped, and within-metric z-scores when metrics are compared on one scale.

## Server scripts

`9x_*.sh` are convenience scripts for syncing to and launching runs on a compute server. Server addresses have been replaced with `<SERVER_IP>` placeholders.
