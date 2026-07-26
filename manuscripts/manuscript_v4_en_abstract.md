# Detection-corrected community dynamics of Chinese birds, 2000–2024

## Structured abstract (≈ 250 words)

**Background.** Citizen-science occurrence records have become a backbone for documenting biodiversity change, yet uneven sampling effort can mistakenly inflate apparent gains in community richness or distribution. Quantifying community-level dynamics for Chinese birds across two-and-a-half decades demands a framework that simultaneously corrects for imperfect detection, accounts for spatial autocorrelation, and propagates uncertainty into derived diversity metrics.

**Methods.** We integrated the China Birdwatching Records Center and curated GBIF/eBird China data (2000–2024), yielding **{n_records}** deduplicated detection events under a source-aware hash. On a 100-km equal-area grid (n = **{n_sites}**) divided into five 5-year primary periods, we fitted a **spatial multi-species dynamic occupancy model (`stMsPGOcc`)** with NNGP spatial random effects (exponential covariance, 5 neighbours) and AR(1) temporal correlation for **{n_species}** species. Posterior ψ samples were propagated jointly into taxonomic (corrected richness, Shannon), functional (CWM, Rao's Q, FEve, FDiv, trait volume) and **probability-weighted Faith's PD** metrics, and into Baselga's β-diversity partitioning (turnover vs nestedness). Drivers were modelled with `brms` using a horseshoe sparsity prior and a spatial Gaussian process (k ∈ {10, 20, 50}). Trait-trend regressions were compared with and without phylogenetic random effects via `loo`.

**Results.** Corrected richness {richness_str}. Adjacent-period β-diversity was {turnover_str}, dominated by turnover rather than nestedness. Of the {n_species} species, {cls_str}. Naive vs occupancy-corrected richness trends flipped direction in {nvc_str} of grids. High-HWI, broad-habitat, low-diet-specialization species disproportionately expanded.

**Conclusions.** Imperfect-detection correction substantially restructures conclusions drawn from raw citizen-science data, even after extensive cleaning. Occupancy-based frameworks should be the default for multi-decadal community-scale inference from such records.

**Keywords.** dynamic occupancy; spatial multi-species model; detection bias; biotic homogenization; Baselga decomposition; probability-weighted PD; citizen science; China.
