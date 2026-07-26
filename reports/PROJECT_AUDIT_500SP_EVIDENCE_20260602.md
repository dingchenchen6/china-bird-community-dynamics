# 中國鳥類群落動態研究本地審查報告

> Date: 2026-06-02  
> Project root: `/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis`  
> Scope: local code, v3/v4 analysis outputs, result tables, figures, with special attention to the 500-species analysis.  
> ARS stage: academic-pipeline intake + integrity-oriented evidence audit.

## Executive Verdict

This project has a strong high-impact core: 25 years of China-scale citizen-science bird records, dynamic occupancy correction, posterior propagation into multiple diversity dimensions, and a 500-species extension. The empirical signal is compelling and internally coherent: occupancy-corrected richness and Shannon diversity rise through time, species-level trends are expansion dominated, beta-diversity shifts from early turnover toward late nestedness, and functional trait space shows mild contraction despite taxonomic gains.

However, the current local repository is not yet submission-ready. The strongest publication framing should treat the **200-species spatial stMsPGOcc model** as the spatially explicit main analysis and the **500-species temporal tMsPGOcc model** as a breadth/robustness extension. The 500-species analysis is explicitly non-spatial in `code_v3/04_run_tMsPGOcc_500sp.R`, and both local model summary tables report `n_chains = 1`. The local `data/derived_v3/` directory does not contain the full model fit objects, so convergence and posterior reproducibility cannot be re-audited locally without resyncing or rerunning.

## What Is Solid Enough To Write

### 1. Data and scope

- Main v3/v4 design: 2000-2024, 5 primary periods, 100 km grid, 1,247 analysed grid cells.
- 200-species spatial result: `table_model_summary_v3_full_200sp_ar1_spatial.csv` reports 200 species, 1,247 sites, 5 periods.
- 500-species temporal result: `table_model_summary_v3_full_500sp_ar1_temporal.csv` reports 500 species, 1,247 sites, 5 periods.
- The local README and scripts document a clear pipeline from source merging to dynamic occupancy, diversity propagation, figures, manuscript, and proposal.

### 2. 200-species spatial model results

The 200-species spatial analysis is the best current basis for spatial maps and mechanistic spatial interpretation.

- Corrected richness increases from 79.15 to 100.83 species per grid on average from P1 to P5.
- Shannon increases from 4.74 to 4.88.
- Inverse Simpson increases from 106.37 to 126.64.
- Functional trait volume declines slightly from 1.383 to 1.365; Rao's Q declines from 1.705 to 1.680; FEve and FDiv rise slightly.
- McTavish phylogenetic diversity increases from 359,581 to 458,111, while MPD increases from 122.28 to 126.30.
- Species trends: 136 expanding, 62 stable, 2 contracting.
- Naive-vs-corrected comparison: 37 of 1,247 grids flip trend direction; median absolute trend difference is 5.994 species per period; mean difference is -8.220.
- Baselga turnover proportions: P1-P2 65.9%, P2-P3 42.8%, P3-P4 53.8%, P4-P5 11.3%.

### 3. 500-species temporal extension

The 500-species analysis strongly supports taxonomic breadth and species-level generality, but because it is temporal and non-spatial it should not be described as a spatial NNGP result.

- Corrected richness increases from 155.04 to 195.55 species per grid on average from P1 to P5.
- Shannon increases from 5.507 to 5.645.
- Inverse Simpson increases from 218.69 to 260.33.
- Functional trait volume declines from 1.581 to 1.565; Rao's Q declines from 1.732 to 1.711; FEve increases very slightly from 0.0690 to 0.0708.
- Trend slopes across grids: corrected richness mean slope 9.563 species per 5-year period; Shannon mean slope 0.0331; inverse Simpson mean slope 9.959.
- Species trends: 258 expanding, 227 stable, 15 contracting.
- Strong expanding species include *Fulica atra*, *Pernis ptilorhynchus*, *Picoides canicapillus*, *Vanellus cinereus*, and *Chlidonias hybrida*.
- Strong contracting/declining species include *Columba livia*, *Seicercus valentini*, *Phylloscopus ogilviegranti*, *Coloeus dauuricus*, and *Lanius tigrinus*.
- Naive-vs-corrected comparison: 25 of 1,247 grids flip trend direction; median absolute trend difference is 8.663 species per period; mean trend difference is -12.037.
- Baselga turnover proportions: P1-P2 75.7%, P2-P3 55.1%, P3-P4 65.3%, P4-P5 19.5%.

### 4. Trait and driver results

Trait/driver tables labelled 500sp contain 500 rows in species-level files, but complete trait/environment covariates are available for only 200 species in the current table joins. Therefore, write these as complete-case mechanism analyses rather than full 500-species mechanism estimates.

Recomputed from `table_species_trend_traits_v3_full_500sp_ar1_temporal_extended.csv`:

- Diet specialization vs trend: Spearman rho = -0.0798, n = 200.
- Habitat breadth vs trend: Spearman rho = 0.3793, n = 200.
- Migration score vs trend: Spearman rho = 0.0802, n = 200.

Recomputed from `table_species_env_trend_v3_full_500sp_ar1_temporal_extended.csv`:

- Cropland cover: rho = 0.3852, n = 200.
- Elevation: rho = 0.3834, n = 200.
- BIO7: rho = 0.3764, n = 200.
- BIO13: rho = 0.3752, n = 200.
- BIO4: rho = 0.3581, n = 200.
- HFI: rho = 0.3576, n = 200.

Variance partitioning for 500sp corrected richness trend:

- Baseline spatial/environment controls: adjusted R2 = 0.1425.
- Climate change: adjusted R2 = 0.0391.
- Land-use change: adjusted R2 = 0.0131.
- Human-pressure change: adjusted R2 = 0.0053.
- Residual: adjusted R2 = 0.6603.

Random forest importance, latest v3 run:

- Top variables: longitude, baseline mean temperature, extreme temperature change, elevation, latitude, impervious surface change.
- Group importance: spatial baseline > land-use change approximately climate change > human-pressure change.

## Code-Level Risks

### P0/P1 risks before submission

1. **500 species is temporal, not spatial.**  
   `code_v3/04_run_tMsPGOcc_500sp.R` lines 8-10 explicitly state the strategy is tMsPGOcc, non-spatial, to avoid spatial NNGP memory failure. Do not write "500-species spatial dynamic occupancy model" unless a new 500sp spatial run is completed.

2. **Local model summaries report one chain.**  
   Both `table_model_summary_v3_full_200sp_ar1_spatial.csv` and `table_model_summary_v3_full_500sp_ar1_temporal.csv` report `n_chains = 1`. The scripts support chain-parallel operation, but the local summary does not prove 4-chain convergence.

3. **Full model fit objects are absent locally.**  
   `data/derived_v3/` contains only small derived RDS files and no full `stMsPGOcc_fit_*` or `tMsPGOcc_fit_*` RDS objects. This blocks local re-audit of posterior arrays, R-hat, ESS, WAIC, PPC, and model metadata.

4. **Spatial phi prior likely needs scale validation.**  
   `code_v3/04_run_stMsPGOcc_main.R` initializes `phi = 3 / mean(pairwise distance km)` but uses `phi.unif = c(0.1, 10)`. At 100 km national scale, the lower bound may be too high unless spOccupancy expects a different distance scale. This should be checked and rerun/sensitivity-tested.

5. **Source-specific detection heterogeneity is not modelled.**  
   Detection model currently uses `log_events + log_duration + has_duration`, but not data source. Because eBird/GBIF and China Birdwatching Records likely differ in protocol and observer structure, add a source or source-effort detection term if the array construction supports it.

6. **Species traits and environmental mechanism analyses are complete-case subsets.**  
   The 500sp trait/environment files have 500 rows, but complete trait/environment covariates for the relevant fields are only 200 species. Mechanism claims must use the complete-case n.

### P2/P3 risks and figure issues

1. **PD display failure in main maps.**  
   The current v3 multidiversity map and trend map show the older Faith's PD panel as nearly all grey. Use McTavish PD/MPD outputs or repair the Jetz tree matching before any top-journal submission.

2. **Temporal synchrony is at boundary 1.0.**  
   Both 200sp and 500sp temporal dynamics summaries show synchrony exactly 1.0 in inspected rows. This is likely a metric implementation/scale issue and should not be emphasized.

3. **Figure aesthetics are not yet top-journal.**  
   Current maps are information-rich but label-dense; axis ticks overlap in the multidiversity timeslice map; RF variable labels are raw variable names; the 500sp analysis has only one exported RF figure.

4. **v4 is a code-review branch, not completed empirical output.**  
   `results_v4/` currently contains audit and draft files but no v4 CSVs, figures, or logs. v4 should be treated as the repair plan unless rerun.

## Recommended Manuscript Architecture

### Best current scientific claim

China's bird communities show a striking paradox: detection-corrected taxonomic diversity and occupancy breadth have increased, especially among generalist and broad-habitat species, while functional trait volume and Rao's Q decline slightly and late-period beta diversity shifts toward nestedness. This points to **expansion without functional differentiation**, a form of hidden homogenization under citizen-science detectability correction.

### Main analysis hierarchy

1. **Main text, primary inference:** 200sp spatial stMsPGOcc results, because these support spatial maps and spatial driver interpretation.
2. **Main text, generality paragraph or SI:** 500sp temporal extension, because it shows the same expansion/homogenization pattern across a much broader species pool.
3. **Mechanism:** complete-case trait/environment analysis, with n = 200 unless the trait joins are expanded.
4. **Sensitivity and validation:** v4 rerun plan, multi-chain diagnostics, source detection term, phi prior/range sensitivity, breeding-season vs year-round, 50 km grid pilot.

## Recommended Figure Set

### Main figures

1. Conceptual workflow: citizen-science records -> dynamic occupancy -> posterior diversity propagation -> beta decomposition -> drivers/traits.
2. 200sp multidiversity time-slice maps, repaired with McTavish PD and cleaner labels.
3. 200sp trend maps for richness, Shannon, trait volume, McTavish PD/MPD.
4. Baselga turnover/nestedness proportion and map.
5. Species winner-loser forest plot with 500sp supplement or split panels by ecological group.
6. Naive vs occupancy-corrected trend comparison.
7. Driver/trait synthesis: RF group importance plus trait effect plot.

### Supplement

- MCMC diagnostics and PPC.
- Breeding-season filtering audit.
- Source-bias maps.
- 3-year vs 5-year and 50 km vs 100 km sensitivity.
- Full 500sp trend table and species list.

## Submission Readiness Gate

Before this can honestly be submitted to Nature Ecology & Evolution, Nature Communications, Global Change Biology, or Ecology Letters:

- Sync or regenerate full model fit objects locally.
- Produce 4-chain R-hat/ESS diagnostics or explicitly label current results as single-chain exploratory.
- Fix or validate the spatial phi prior/range scale.
- Add source-specific detection heterogeneity or justify why it is not identifiable.
- Repair PD maps and switch to verified McTavish phylogenetic metrics.
- Generate a 500sp figure set beyond RF importance.
- Freeze a reproducible v4 run with logs, sessionInfo, code hash, and data manifest.

## Verified External Anchors Used For Rewriting

- spOccupancy package paper: Doser et al., Methods in Ecology and Evolution, DOI `10.1111/2041-210X.13897`.
- Baselga beta-diversity partitioning: Baselga 2010, Global Ecology and Biogeography, DOI `10.1111/j.1466-8238.2009.00490.x`.
- Closure assumption in occupancy: Rota et al. 2009, Journal of Applied Ecology, DOI `10.1111/j.1365-2664.2009.01734.x`.
- AVONET: Tobias et al. 2022, Ecology Letters, DOI `10.1111/ele.13898`.
- EltonTraits: Wilman et al. 2014, Ecology, DOI `10.1890/13-1917.1`.
- Avian diet specialization: Morelli et al. 2021, Conservation Letters, DOI `10.1111/conl.12795`.
- BirdTree/Jetz phylogeny: Jetz et al. 2012, Nature, DOI `10.1038/nature11631`.
- China 30 m land cover dataset: Yang & Huang 2021, Earth System Science Data, DOI `10.5194/essd-13-3907-2021`.
- Annual terrestrial Human Footprint: Mu et al. 2022, Scientific Data, DOI `10.1038/s41597-022-01284-8`.
