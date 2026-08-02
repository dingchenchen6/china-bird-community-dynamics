# Detection-corrected dynamics of Chinese bird communities, 2000–2024

Code, result tables, figures and manuscripts for a national-scale analysis of Chinese bird community change, in which citizen-science records are corrected for imperfect detection with spatial multispecies dynamic occupancy models and the occupancy posterior is propagated into taxonomic, functional, phylogenetic and beta-diversity space.

**Headline result.** Detection-corrected richness rose by roughly a quarter (79.2 → 100.8 species per 100-km grid, 2000–2024), while functional trait volume and Rao's quadratic entropy did not rise, temporal beta diversity shifted from turnover to nestedness, and expanding species were overwhelmingly broad-habitat generalists — taxonomic expansion without commensurate functional differentiation.

> ⚠️ **Work in progress.** The associated manuscript is unpublished and under preparation. Results in `results/` and `manuscripts/` may change. Please do not cite without contacting the author.

## Study design

| Item | Value |
|---|---|
| Period | 2000–2024, five 5-year periods |
| Grid | 100 km equal-area, 1,247 cells with data |
| Primary model | `spOccupancy::stMsPGOcc` — 200 species, spatial NNGP + AR(1) |
| Breadth extension | `spOccupancy::tMsPGOcc` — 500 species, temporal (non-spatial), 10 latent factors |
| MCMC | 4 chains, 5,000 burn-in, thin 2, 1,250 post-burn-in draws per chain |
| Data sources | China Birdwatching Records Center + curated GBIF / eBird China records |
| Traits & phylogeny | AVONET, EltonTraits, IUCN-derived habitat breadth, global bird tree |

## Where to start

| If you want to… | Read |
|---|---|
| Re-run the analysis | [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) — environment, inputs, pipeline order, runtimes |
| Run the mechanism & conservation analyses on a compute server | [`server_run/README_KIMI_SERVER_RUN.md`](server_run/README_KIMI_SERVER_RUN.md) — self-contained, one-command driver |
| Understand the result tables | [`results/README.md`](results/README.md) — full data dictionary |
| Understand the code | [`code/README.md`](code/README.md) — stage map and conventions |
| Judge how solid the study is | [`reports/RESEARCH_SYSTEMATIC_AUDIT_20260602.md`](reports/RESEARCH_SYSTEMATIC_AUDIT_20260602.md) — end-to-end audit |
| Read the paper | [`manuscripts/README.md`](manuscripts/README.md) |

## Repository structure

```
code/         R pipeline                        → code/README.md
  code_v3/    current main pipeline (data prep → occupancy → posterior propagation → figures)
  code_v4/    audit / repair branch
  code_v2/    earlier version, retained for comparison
results/      result tables (CSV)               → results/README.md (data dictionary)
  results_v3/ current main results (4-chain)
  results_v4/ audit-stage outputs
manuscripts/  manuscript drafts, cover letters, abstracts (Markdown + DOCX)
reports/      methodological audits, simulated peer review, rerun guides
figures/      current main-line figures (PNG/PDF)
server_run/   one-command driver + self-contained server instructions
REPRODUCIBILITY.md   how to reproduce, end to end
```

### Mechanism and conservation analyses (scripts 30–32)

Three analyses address why richness rose and what it means for conservation:

| Script | Question | Key output |
|---|---|---|
| `30_range_expansion_mechanism.R` | Is rising richness range-edge expansion or in-range infilling? Climate-driven? | Range decomposition, centroid shifts, community temperature index, guild attribution |
| `31_protected_area_effectiveness.R` | Do China's nature reserves slow functional homogenization? | Propensity-score-matched ATT, difference-in-differences, representation gaps |
| `32_conservation_prioritization.R` | Does richness-led planning miss the areas that prevent homogenization? | Three prioritization scenarios, spatial mismatch, protected-area shortfall |

Protected-area boundaries come from China Nature Reserve Specimen Resource Sharing Platform (2024), *List and Vector Boundaries of Nature Reserves in China*, Zenodo, https://doi.org/10.5281/zenodo.14875797 (CC-BY-4.0); the data are not redistributed here.

### Pipeline entry points (`code/code_v3/`)

| Stage | Script |
|---|---|
| Merge & deduplicate records | `01_merge_birdwatch_ebird.R` |
| Build survey histories | `02_build_survey_history.R` |
| Environmental covariates | `03_prepare_environment.R` (+ `03b`–`03e`) |
| Fit occupancy models | `04_run_stMsPGOcc_main.R`, `04_run_tMsPGOcc_500sp.R` |
| Posterior diversity propagation | `05_postprocess_diversity_extended.R` |
| Diagnostics & PPC | `05b_mcmc_diagnostic_plots.R`, `05c_ppc_bayesian_pvalue.R` |
| Functional trend posteriors | `05f_functional_trend_pdecline.R` |
| Figures | `06_figures_publication.R`, `21`–`24_*` |
| Sensitivity analyses | `15`–`17_sensitivity_*.R` |

## Data availability

Raw occurrence records are **not** included here and are subject to source-provider terms (China Birdwatching Records Center, GBIF, eBird). Large intermediate objects (posterior samples, model fits; tens of GB) are also excluded — see `.gitignore`. This repository contains the code, derived result tables, figures and manuscripts needed to follow and audit the analysis. Derived data will be archived with a DOI on publication.

## Reproducing

The pipeline is written in R and depends principally on `spOccupancy`, `brms`, `cmdstanr`, `sf`, `terra`, `vegan`, `ape` and `ggplot2`. Scripts read shared configuration from `code/code_v3/00_config.R` and resolve paths through `utils_paths.R`; set `V3_CODE_DIR` to point at the code directory. Large model fits must be regenerated (stage 04) before the postprocessing stages will run.

## Citation & contact

Chenchen Ding, Peking University — jialeding1220@gmail.com
Website: https://dingchenchen6.github.io/chenchen-ding/ · GitHub: https://github.com/dingchenchen6

## License

Code is released under the MIT License (`LICENSE`). Result tables and figures are © the author; please contact before reuse while the manuscript is unpublished.
