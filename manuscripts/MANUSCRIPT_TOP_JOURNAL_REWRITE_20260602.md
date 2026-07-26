# Expansion Without Differentiation: Occupancy-Corrected Citizen Science Reveals Functional Homogenization of Chinese Bird Communities, 2000-2024

**Short title:** Detection-corrected bird community change in China  
**Article type:** Research Article / Nature Ecology & Evolution, Nature Communications, Global Change Biology, or Ecology Letters target  
**Evidence status:** This draft is based on local v3 result tables and figures audited on 2026-06-02. It is written as a high-impact manuscript draft, but the submission version must replace single-chain model summaries with verified multi-chain convergence diagnostics and sync or regenerate the full model fit objects.

## Abstract

Citizen-science biodiversity records now provide continental-scale evidence for ecological change, but uneven observation effort can turn expanded sampling into apparent biological gain. Here we integrated China Birdwatching Records Center observations with curated GBIF/eBird China records from 2000 to 2024 and analysed bird community dynamics on a 100-km grid across five 5-year periods. We fitted a spatial multispecies dynamic occupancy model for 200 species, correcting imperfect detection and propagating posterior occupancy probabilities into taxonomic, functional and phylogenetic diversity, beta-diversity partitioning, and species-level trends. We further used a 500-species temporal dynamic occupancy extension to test whether the principal patterns generalize across a broader species pool.

Occupancy-corrected richness increased strongly in the 200-species spatial model, from 79.15 to 100.83 species per grid on average, and the 500-species temporal extension showed a parallel increase from 155.04 to 195.55. Shannon and inverse Simpson diversity also increased, yet functional trait volume and Rao's Q declined slightly, indicating taxonomic expansion without functional differentiation. Species-level trends were expansion dominated: 136 of 200 spatial-model species and 258 of 500 temporal-extension species expanded, whereas only 2 and 15 species, respectively, contracted. Temporal beta diversity shifted from turnover-dominated early intervals to strong nestedness in 2020-2024, with turnover proportions falling from 65.9% to 11.3% in the 200-species model and from 75.7% to 19.5% in the 500-species extension. Habitat breadth was positively associated with species occupancy trends, whereas diet specialization was weakly negative in the complete-case trait subset.

Our results reveal a paradox of contemporary avian change in China: more species are expected per grid after detection correction, but communities increasingly favour broad-habitat, functionally redundant assemblages. Occupancy correction changed trend magnitudes and flipped trend directions in 2-3% of grids, showing that citizen-science biodiversity inference requires explicit modelling of detectability before ecological interpretation.

**Keywords:** dynamic occupancy, citizen science, birds, China, detection bias, biotic homogenization, functional diversity, beta diversity, Baselga partitioning, AVONET

## Introduction

The central problem in biodiversity change research is no longer only whether enough observations exist. It is whether the observations can distinguish ecological change from changing observation itself. Citizen-science bird data are uniquely powerful because they provide repeated, spatially extensive records across decades, but they are also strongly non-random: observers visit accessible places, effort changes over time, platforms differ in protocol, and detection varies among species. For countries such as China, where birdwatching participation and digital reporting have expanded rapidly, raw species counts can easily confuse increased detectability with community recovery.

Dynamic occupancy models offer a principled response to this problem. By separating latent occurrence from detection probability, they allow community inference from detection-nondetection data while accounting for imperfect observation. Spatial multispecies implementations extend this logic further, sharing information among species, borrowing strength across space, and propagating uncertainty into derived biodiversity metrics. Yet few studies have used such models to connect nationwide citizen-science records with multidimensional community change, especially in regions undergoing rapid climate warming, land-use transformation, urban expansion, and conservation investment.

China is an unusually important test case. It contains major East Asian flyways, high mountain systems, arid western basins, humid subtropical forests, intensively transformed agricultural lowlands, and densely urbanized coastal regions. These gradients create a strong expectation of spatially heterogeneous bird community change. At the same time, observer coverage is uneven and historically shifting. A national analysis therefore requires a framework that can correct detection, represent spatial structure, and ask whether apparent increases in richness correspond to genuine ecological differentiation or to homogenized gain by widespread generalists.

Here we analyse Chinese bird community dynamics from 2000 to 2024 using a spatial multispecies dynamic occupancy framework and posterior propagation into taxonomic, functional, phylogenetic, and beta-diversity metrics. We ask five questions. First, how have occupancy-corrected taxonomic and multidimensional diversity changed across China? Second, are temporal changes driven by species replacement or nested subsets of communities? Third, which species are expanding or contracting after detection correction? Fourth, do functional traits and broad environmental gradients explain species-level trends? Fifth, how much do conclusions change relative to naive occurrence-based estimates? We use a 200-species spatial model as the primary spatial inference and a 500-species temporal extension as a breadth test of the main biological signal.

## Results

### Modelled scope and evidence hierarchy

The primary spatial analysis modelled 200 bird species across 1,247 100-km grid cells and five 5-year periods from 2000-2004 to 2020-2024. The extension analysis modelled 500 species across the same grid-period structure using a temporal multispecies dynamic occupancy model without spatial NNGP terms. We therefore interpret the 200-species model as the basis for spatial pattern and driver inference, and the 500-species model as a broader test of taxonomic generality.

### Taxonomic diversity increased after occupancy correction

Occupancy-corrected richness increased monotonically in the 200-species spatial model. Mean grid-level corrected richness rose from 79.15 species in P1 to 82.23 in P2, 87.02 in P3, 90.53 in P4, and 100.83 in P5. Median corrected richness followed the same pattern, increasing from 80.35 to 102.74 species per grid. Shannon diversity increased from 4.74 to 4.88 and inverse Simpson diversity increased from 106.37 to 126.64.

The 500-species temporal extension showed the same direction at larger magnitude. Mean corrected richness increased from 155.04 to 195.55 species per grid, median corrected richness from 149.08 to 191.24, Shannon from 5.507 to 5.645, and inverse Simpson from 218.69 to 260.33. Across grids, the mean corrected-richness slope was 9.56 species per 5-year period in the 500-species extension, with a median grid slope of 8.87.

### Functional diversity did not keep pace with taxonomic gain

Taxonomic increases were not accompanied by expansion of functional trait space. In the 200-species spatial model, functional trait volume declined from 1.383 to 1.365, while Rao's Q declined from 1.705 to 1.680. FEve and FDiv increased slightly, from 0.164 to 0.171 and from 0.691 to 0.699, respectively. This combination suggests more even filling of an existing trait space rather than expansion into new functional extremes.

The 500-species temporal extension sharpened the same interpretation. Trait volume declined from 1.581 to 1.565 and Rao's Q from 1.732 to 1.711, whereas FEve increased only slightly from 0.0690 to 0.0708. These changes point to a national-scale pattern of taxonomic expansion without functional diversification.

### Phylogenetic metrics increased where tree matching succeeded

McTavish phylogenetic outputs in the 200-species spatial model showed increasing phylogenetic diversity and mean pairwise phylogenetic distance. Probability-weighted phylogenetic diversity increased from 359,581 to 458,111, and MPD increased from 122.28 to 126.30. Earlier Jetz-tree PD panels in the local figures appear incomplete or greyed out and should not be used in the submission figure set until tree matching is repaired or replaced by the verified McTavish outputs.

### Species-level trends were expansion dominated

In the 200-species spatial model, 136 species were classified as expanding, 62 as stable, and 2 as contracting. The fastest expanding species included *Fulica atra*, *Pernis ptilorhynchus*, *Picoides canicapillus*, *Chlidonias hybrida*, and *Vanellus cinereus*. The strongest contraction was estimated for *Columba livia*, with *Periparus rubidiventris* also classified as contracting.

The 500-species temporal extension preserved this expansion-dominated structure but revealed more heterogeneity in the broader species pool: 258 species expanded, 227 were stable, and 15 contracted. The fastest expanding species were *Fulica atra*, *Pernis ptilorhynchus*, *Picoides canicapillus*, *Vanellus cinereus*, *Chlidonias hybrida*, *Ardea intermedia*, *Himantopus himantopus*, *Corvus dauuricus*, *Luscinia svecica*, and *Sturnus vulgaris*. The most negative trends were estimated for *Columba livia*, *Seicercus valentini*, *Phylloscopus ogilviegranti*, *Coloeus dauuricus*, *Lanius tigrinus*, *Certhia hodgsoni*, *Locustella thoracica*, and *Lophophanes dichrous*.

### Temporal beta diversity shifted toward nestedness in the latest interval

Probability-weighted Baselga partitioning showed strong temporal variation in the relative importance of turnover and nestedness. In the 200-species spatial model, turnover accounted for 65.9% of beta diversity from P1 to P2, 42.8% from P2 to P3, 53.8% from P3 to P4, and only 11.3% from P4 to P5. In the 500-species temporal extension, turnover accounted for 75.7%, 55.1%, 65.3%, and 19.5% over the same intervals.

This late decline in turnover proportion indicates that recent changes are less about balanced species replacement and more about communities becoming nested subsets of each other. Combined with increasing corrected richness and mild functional contraction, this pattern is consistent with broad-habitat species accumulating across many grids while distinct local combinations weaken.

### Traits and environmental associations point to broad-habitat expansion

Trait and environmental mechanism analyses were currently complete for 200 species in the 500sp-labelled trait/environment joins. Within this complete-case set, habitat breadth was positively associated with species trend slopes (Spearman rho = 0.379), whereas diet specialization was weakly negative (rho = -0.080). Migration score had little association (rho = 0.080). These results support the interpretation that broad-habitat species are more likely to expand.

Species trend slopes were positively correlated with cropland cover (rho = 0.385), elevation (rho = 0.383), temperature annual range (BIO7, rho = 0.376), wettest-month precipitation (BIO13, rho = 0.375), temperature seasonality (BIO4, rho = 0.358), human footprint (rho = 0.358), and built land cover (rho = 0.339). These are spatial association results, not causal effects of environmental change.

Variance partitioning of the 500-species corrected-richness trend attributed the largest pure component to baseline spatial/environmental structure (adjusted R2 = 0.1425), followed by climate change (0.0391), land-use change (0.0131), and human-pressure change (0.0053), with substantial residual variance (0.6603). Random forest permutation importance similarly ranked spatial baseline variables highest, especially longitude, followed by baseline winter temperature, extreme temperature change, elevation, latitude, and impervious surface change.

### Occupancy correction changed naive trend inference

Naive occurrence trends and occupancy-corrected trends were directionally similar for most grids, but correction changed trend magnitude substantially. In the 200-species spatial model, 37 of 1,247 grids (3.0%) flipped trend direction after correction, with a median absolute trend difference of 5.99 species per period and a mean difference of -8.22. In the 500-species temporal extension, 25 grids (2.0%) flipped direction, with a median absolute trend difference of 8.66 species per period and a mean difference of -12.04.

Thus, detectability correction did not simply rescale raw richness: it changed the magnitude of inferred change across many grids and reversed interpretation in a small but non-negligible set of locations.

## Discussion

### A national signal of expansion without differentiation

The most striking result is not simply that occupancy-corrected richness increased. It is that this increase occurred while functional trait volume and Rao's Q declined slightly and late-period beta diversity became increasingly nested. These joint results suggest an "expansion without differentiation" pattern: more species are expected per grid, but the species being added do not expand functional trait space in proportion to their taxonomic contribution. Such a pattern is consistent with biotic homogenization, where widespread or broad-habitat species increasingly dominate community change.

This finding matters because raw citizen-science richness trends could easily be read as biodiversity recovery. Occupancy correction makes the result more credible, but the multidimensional analysis makes it more cautious. More species per grid does not necessarily mean more distinctive communities, more functional complementarity, or reduced conservation concern.

### Why the 500-species extension changes the story

The 500-species temporal extension strengthens the manuscript because it shows that the expansion signal is not an artefact of a narrow 200-species spatial subset. In the broader species pool, corrected richness, Shannon diversity, and inverse Simpson diversity all increased, and the expansion-dominated species trend distribution persisted. At the same time, the 500-species model revealed more contracting species and a stronger median correction effect in naive-vs-occupancy comparisons.

The limitation is equally important: the 500-species run is temporal and non-spatial. It should not be used to infer spatial random effects or spatially explicit driver mechanisms. Its role is to show breadth and consistency, while the 200-species spatial model supports maps and spatial hypotheses.

### Turnover-to-nestedness shift as an early warning

The transition from early turnover dominance to late nestedness dominance is a key ecological result. Early intervals involved more species replacement among periods, whereas the P4-P5 interval was overwhelmingly nestedness dominated. This may indicate that recent communities are becoming more similar through shared gains of widespread species or through uneven local retention of a common species pool. Because corrected richness increased, this is not a simple local extinction story. Instead, it suggests convergence in community composition under broad regional changes in observation-corrected occupancy.

### Trait filters identify likely winners

The positive association between habitat breadth and species trends provides a simple biological mechanism for the observed pattern. Species able to occupy more habitat categories are likely to exploit heterogeneous, modified, or newly surveyed landscapes. The weak negative relationship with diet specialization is consistent with specialist vulnerability, although the current complete-case n is only 200 species and should be expanded before strong trait claims are made. The top expanding species include several waterbirds, open-country species, and human-associated or mobile taxa, aligning with a broad-habitat interpretation.

### Environmental associations are spatially structured, not causal proof

Environmental correlations and random forest importance emphasize spatial baseline structure, longitude, winter temperature, elevation, and climate-change variables. These results should be interpreted as spatial associations with occupancy trends rather than causal attribution. Much variation remains unexplained, and spatial structure can absorb unmeasured processes such as observer expansion, regional conservation activity, habitat restoration, agricultural intensification, and platform-specific reporting patterns. The submission version should report marginal and conditional R2 for spatial driver models and avoid causal language unless lagged or quasi-experimental analyses are added.

### Why occupancy correction changes conservation interpretation

The naive-vs-corrected comparison shows that detection correction mainly changes trend magnitude, with a smaller number of direction flips. This is exactly what one would expect when increased sampling effort and true biological change move in the same direction in many regions. Correcting detection does not erase the national richness increase, but it prevents overinterpreting raw increases as ecological recovery. For conservation planning, the most important grids are those where naive and corrected trends disagree, and those where corrected richness increases while functional diversity declines or beta diversity becomes nested.

### Limitations and submission-critical improvements

Several limitations must be resolved or explicitly framed before submission. First, local model summaries report single-chain runs, and full model fit objects are absent from the local derived-data directory. The final paper requires multi-chain R-hat/ESS diagnostics or a transparent explanation of exploratory single-chain status. Second, the 500-species analysis is non-spatial; it is powerful as a breadth extension but not as a spatial mechanism model. Third, source-specific detection heterogeneity should be modelled or tested because eBird/GBIF and China Birdwatching Records differ in protocol. Fourth, the spatial phi prior should be scale-checked for 100-km national distances. Fifth, older PD map panels appear incomplete and should be replaced with verified McTavish phylogenetic metrics.

## Methods

### Data integration

We combined records from the China Birdwatching Records Center and curated GBIF/eBird China records for 2000-2024. Records were cleaned to harmonize species names, dates, geographic coordinates, and source identifiers. The v3 workflow uses a source-aware deduplication strategy based on species, date, rounded coordinates, source, and observer identity. This strategy is designed to avoid collapsing different source platforms while removing within-source duplicate observations. The final submission should include a reproducible deduplication audit and source-wise annual coverage table.

### Spatial and temporal design

We aggregated observations to a 100-km grid and five 5-year primary periods: 2000-2004, 2005-2009, 2010-2014, 2015-2019, and 2020-2024. The primary spatial analysis used 1,247 grid cells with data. Within-period years were used as repeated temporal occasions where supported by the survey-history construction. This is a relaxed-closure design and should be justified with sensitivity analysis using shorter temporal windows.

### Detection and occupancy models

The 200-species primary analysis used a spatial multispecies dynamic occupancy model implemented with `spOccupancy::stMsPGOcc`. The occupancy component included climatic, topographic, habitat, human-footprint, land-cover, spatial-coordinate, and scaled-time covariates. The detection component included `log_events`, `log_duration`, and `has_duration`, where the latter marks missing or unavailable duration information. The model included AR(1) temporal structure and NNGP spatial random effects with exponential covariance.

The 500-species extension used `spOccupancy::tMsPGOcc` with AR(1) temporal structure and 10 latent factors. This model intentionally omitted spatial NNGP structure to avoid memory failure for 500 species. We therefore use it as a temporal breadth extension rather than as a spatial inference engine.

### Posterior diversity propagation

For each posterior occupancy draw, site, and period, we calculated expected species richness as the sum of species occupancy probabilities, along with Shannon diversity, inverse Simpson diversity, functional trait volume, trait dispersion, Rao's Q, FEve, FDiv, and phylogenetic metrics where tree matching succeeded. Functional traits were drawn from AVONET, EltonTraits, IUCN-derived habitat breadth, and imputed trait tables. For phylogenetic diversity, the submission version should use the verified probability-weighted McTavish outputs unless Jetz tree matching is repaired.

### Temporal beta diversity

We computed probability-weighted Baselga beta-diversity partitioning for adjacent periods. For two occupancy probability vectors, expected shared occupancy was defined as the sum of pairwise minima across species. Turnover and nestedness components were then calculated using Baselga-style formulas generalized to probability-weighted inputs. A mathematical supplement should be included in the submission version to prove bounds and decomposition behaviour under probabilistic inputs.

### Species trends and classification

Species-level trends were estimated from posterior occupancy trajectories across periods. Species were classified as expanding, stable, or contracting based on posterior trend direction and credible intervals. Grid-level diversity trends were summarized with OLS and robust trend estimators where available; the final submission should prioritize Theil-Sen slopes and Mann-Kendall tests for 5-period time series.

### Trait and driver analysis

Trait-trend and environmental-trend associations were calculated for complete-case species subsets. Habitat breadth, diet specialization, migration score, climate variables, elevation, human footprint, built cover, and cropland cover were tested against species occupancy trend slopes using Spearman correlations. Community trend drivers were analysed with variance partitioning and random forest permutation importance. Because these analyses are spatially structured and observational, we interpret them as associations rather than causal effects.

### Software and reproducibility

The project is written primarily in R and uses `spOccupancy`, `brms`, `cmdstanr`, `ape`, `vegan`, `sf`, `terra`, `ggplot2`, and related packages. The submission package should include a locked sessionInfo file, code hash, data manifest, model-fit checksums, and a Zenodo or institutional archive DOI.

## Main Figure Plan

**Figure 1. Study design and inference framework.** Record integration, source-aware deduplication, dynamic occupancy correction, posterior diversity propagation, beta partitioning, trait/driver analysis, and naive-vs-corrected comparison.

**Figure 2. Occupancy-corrected multidiversity dynamics.** 200-species spatial maps across five periods for corrected richness, Shannon, functional trait volume, Rao's Q, FEve/FDiv, and verified McTavish phylogenetic metrics. Current PD panels should be repaired before use.

**Figure 3. Trend geography.** Spatial maps of corrected-richness, Shannon, functional trait volume, and phylogenetic trend z-scores with uncertainty masks.

**Figure 4. Temporal beta-diversity reorganization.** Turnover and nestedness proportions across adjacent periods, with the P4-P5 nestedness shift highlighted.

**Figure 5. Species winners and losers.** Posterior trend forest plot for leading expanding and contracting species, with a 500-species supplementary panel.

**Figure 6. Traits and environmental associations.** Habitat breadth and diet specialization effects, plus RF group importance and key variables.

**Figure 7. Naive versus occupancy-corrected trends.** Paired grid-level trend comparison, highlighting direction flips and magnitude differences.

## Data and Code Availability

Code is available in the local repository under `code_v3/` and `code_v4/`. The current local archive lacks the full model fit RDS objects for the 200-species and 500-species analyses; these must be synced or regenerated before public archiving. Raw records are subject to source-provider terms. Derived tables and figures should be archived with a DOI after final verification.

## Acknowledgements

We acknowledge the China Birdwatching Records Center, eBird and the Cornell Lab of Ornithology, GBIF, AVONET contributors, EltonTraits contributors, IUCN Red List, WorldClim, CLCD, Human Footprint dataset authors, and the developers of `spOccupancy` and associated open-source ecological modelling tools.

## References

Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19, 134-143. https://doi.org/10.1111/j.1466-8238.2009.00490.x

Doser, J. W., Finley, A. O., Kery, M., & Zipkin, E. F. (2022). spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. *Methods in Ecology and Evolution*. https://doi.org/10.1111/2041-210X.13897

Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. *Nature*, 491, 444-448. https://doi.org/10.1038/nature11631

Morelli, F., Benedetti, Y., Hanson, J. O., & Fuller, R. A. (2021). Global distribution and conservation of avian diet specialization. *Conservation Letters*, 14, e12795. https://doi.org/10.1111/conl.12795

Mu, H., Li, X., Wen, Y., Huang, J., Du, P., Su, W., Miao, S., & Geng, M. (2022). A global record of annual terrestrial Human Footprint dataset from 2000 to 2018. *Scientific Data*. https://doi.org/10.1038/s41597-022-01284-8

Rota, C. T., Fletcher, R. J. Jr., Dorazio, R. M., & Betts, M. G. (2009). Occupancy estimation and the closure assumption. *Journal of Applied Ecology*, 46, 1173-1181. https://doi.org/10.1111/j.1365-2664.2009.01734.x

Tobias, J. A., Sheard, C., Pigot, A. L., et al. (2022). AVONET: Morphological, ecological and geographical data for all birds. *Ecology Letters*, 25, 581-597. https://doi.org/10.1111/ele.13898

Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M. M., & Jetz, W. (2014). EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. *Ecology*, 95, 2027. https://doi.org/10.1890/13-1917.1

Yang, J., & Huang, X. (2021). The 30 m annual land cover dataset and its dynamics in China from 1990 to 2019. *Earth System Science Data*, 13, 3907-3925. https://doi.org/10.5194/essd-13-3907-2021
