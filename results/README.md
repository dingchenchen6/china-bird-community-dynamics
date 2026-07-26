# Result tables — data dictionary

All tables are CSV. `results_v3/` holds the current main results (4-chain models); `results_v4/` holds audit-stage outputs.

## Naming convention

```
table_<content>_<run_label>[_extended].csv
```

`<run_label>` identifies the model run and is one of:

| Run label | Model | Species | Spatial? | Role |
|---|---|---|---|---|
| `v3_full_200sp_ar1_spatial` | `stMsPGOcc` | 200 | yes (NNGP) | **primary inference** — all spatial maps and mechanisms |
| `v3_full_500sp_ar1_temporal` | `tMsPGOcc` | 500 | no | breadth / generality test only |
| `v3_pilot_60sp_ar1_spatial` | `stMsPGOcc` | 60 | yes | pilot run |

Suffixes: `_extended` = extended postprocessing (adds functional/phylogenetic metrics); `_4chain` = 4-chain diagnostics; `_200draws` = 200-posterior-draw variant.

⚠️ Never use a `*_500sp_ar1_temporal` table for a spatial claim — that model has no spatial random effects.

## Shared column patterns

Most posterior-summary tables share this shape:

| Column | Meaning |
|---|---|
| `mean`, `sd` | posterior mean and standard deviation |
| `q025`, `median`, `q975` | 2.5%, 50%, 97.5% posterior quantiles (95% credible interval) |
| `grid_cell` | 100-km grid cell identifier |
| `period` | `P1`–`P5` = 2000–2004, 2005–2009, 2010–2014, 2015–2019, 2020–2024 |
| `period_pair` | adjacent-period transition, e.g. `P4_P5` |
| `metric` | diversity metric name (see below) |

**Metric vocabulary:** `corrected_richness` (detection-corrected species richness), `shannon`, `inv_simpson`, `trait_volume` (functional richness / trait-space volume), `trait_dispersion`, `rao_q` (Rao's quadratic entropy), `feve`/`feve_fund` (functional evenness), `fdiv`/`fdis_prob` (functional divergence/dispersion), `fric_prob` (probabilistic functional richness), `pd_prob` / `pd_prob_mctavish` (probability-weighted phylogenetic diversity), `mpd_prob` / `mpd_prob_mctavish` (mean pairwise phylogenetic distance).

## Key tables

### Model & convergence
| Table | Columns | Notes |
|---|---|---|
| `table_model_summary_*` | run metadata | species, sites, periods, chain settings |
| `table_convergence_diagnostics_*_4chain` | `parameter`, `group`, `rhat`, `ess` | 38 community-level parameters; max R̂ = 1.039, min ESS = 70 |

### Diversity (per grid × period)
| Table | Columns |
|---|---|
| `table_diversity_summary_*` | `mean`, `sd`, `q025`, `median`, `q975`, `metric`, `period`, `grid_cell` |
| `table_community_metrics_with_cri_*` | `grid_cell`, `period`, `metric`, `value_mean`, `value_l95`, `value_u95`, `block_label` |
| `table_diversity_wide_*` | wide form, one column per metric |

### Trends (per grid)
| Table | Columns |
|---|---|
| `table_trend_summary_*` | `mean`, `sd`, `q025`, `median`, `q975`, `metric`, `method`, `grid_cell` — `method` distinguishes OLS vs Theil–Sen |
| `table_mann_kendall_*` | Mann–Kendall trend tests |
| `table_naive_vs_corrected_*` | naive occurrence trend vs occupancy-corrected trend, incl. direction flips |

### Beta diversity / homogenization
| Table | Columns |
|---|---|
| `table_baselga_summary_*` | `mean`, `sd`, `q025`, `median`, `q975`, `metric` (`beta_sor`/`beta_sim`/`beta_sne`/`prop_turnover`), `period_pair`, `grid_cell` |
| `table_baselga_global_*` | `period_pair`, `prop_turnover_mean`, `prop_turnover_q025`, `prop_turnover_q975`, `n_grids` — national turnover proportion |
| `table_homogenization_trend_*` | `period`, `mean_sorensen`, `sd_sorensen` — among-grid spatial beta ⚠️ currently `NA`, see `reports/RESEARCH_SYSTEMATIC_AUDIT_20260602.md` |

### Species level
| Table | Columns |
|---|---|
| `table_species_trend_traits_*` | `species`, `trend_class`, `trend_grade`, `trend_slope`, `p_positive`, `p_negative`, `mk_tau`, `mk_p`, `diet_specialization`, `habitat_breadth`, `Habitat`, `Habitat.Density`, `Migration`, `Trophic.Level`, `Primary.Lifestyle`, `migration_score` |
| `table_species_trend_classify_*` | expanding / stable / contracting classification |
| `table_species_winners_losers_*` | leading expanding and contracting species |
| `table_species_env_trend_*` | species trend vs environmental covariates |

### Drivers
| Table | Columns |
|---|---|
| `table_varpart_<metric>_*` | `fraction`, `Df`, `R.square`, `adj_r2`, `Testable` — variation partitioning |
| `table_rf_importance_group_summary_trend` | `group`, `mean`, `median`, `sd`, `q_lo`, `q_hi`, `label` — random-forest group importance |
| `table_rf_importance_variable_summary_trend` | per-variable permutation importance |
| `table_brms_driver_trend_change_coefficients_<metric>_*` | brms driver-regression coefficients |
| `table_env_trend_correlation_*`, `table_trait_trend_correlation_*` | Spearman correlations |

### Validation
| Table | Columns |
|---|---|
| `table_effort_confound_controls_*` | `control`, `parameter`, `mean`, `sd`, `q025`, `median`, `q975`, `n_grids` — detection–effort slope, effort-saturated subset, fixed-effort counterfactual, source term |
| `table_breeding_filter_audit` | breeding-season sensitivity |
| `table_vif_screening_*` | collinearity screening |
| `table_stability_summary_*` | temporal stability metrics |

### Inputs / audit
`table_01_merge_summary`, `table_02_survey_summary*` (record merging and survey-history construction), `table_grid_environment*`, `table_climate_change`, `table_landuse_change`, `table_hfi_change` (covariates), `table_traits_extended` (trait table).

## Caveats

- Tables labelled `500sp` for traits/environment contain 500 rows but **complete covariates exist for 200 species only** — mechanism analyses are complete-case (n = 200).
- `table_temporal_dynamics_summary_*` reports synchrony = 1.0, a known metric-implementation artefact; do not interpret.
