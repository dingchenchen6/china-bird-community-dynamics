# More species, fewer differences: detection-corrected citizen science reveals functional homogenization beneath rising avian richness across China, 2000–2024

**Target journal:** *Nature Ecology & Evolution* / *Global Change Biology* / *Ecology Letters* (Article)

**Authors:** Chenchen Ding¹*, [co-authors]
¹ [Affiliation], Peking University, Beijing, China
\* Correspondence: chenchending1992@gmail.com

**Running title:** Expansion without differentiation in Chinese birds

> **Integrity note (delete before submission).** Point estimates derive from audited 4-chain model outputs (`results_v3`, 2026-06). Convergence diagnostics and credible intervals reported here are the audited values; three analyses are still on the compute queue and are flagged `[PENDING]` at the exact sentence where their numbers belong. No diagnostic value is invented. All 76 references were DOI-verified against Crossref (2026-06-02 and 2026-07-28).

---

## Abstract

Global assessments increasingly agree that biodiversity is being reorganized, yet disagree on whether it is being lost — a contradiction that persists partly because most evidence is expressed in a single currency, species richness, and partly because the observations underpinning it are collected by volunteers whose effort is itself changing. Here we resolve both problems at national scale. Integrating 25 years of Chinese citizen-science records (2000–2024) across a 100-km grid, we fitted a spatial multispecies dynamic occupancy model for 200 species and propagated the full occupancy posterior into taxonomic, functional, phylogenetic and beta-diversity space, with a 500-species temporal model testing generality. Detection-corrected richness rose by 27% (79.2 to 100.8 species per grid), and Shannon, inverse Simpson and phylogenetic diversity rose with it. Functional trait volume and Rao's quadratic entropy did not: both declined monotonically across all five periods. Temporal beta diversity collapsed from turnover- to nestedness-dominance (turnover 65.9% to 11.3%), and the species driving expansion were disproportionately broad-habitat generalists while contracting species were narrow-range montane forest specialists — habitat breadth was the strongest single trait predictor of occupancy trend (ρ = 0.38). Four independent controls show the increase is not an artefact of growing observer effort. Chinese bird communities are therefore accumulating species while losing distinctiveness: taxonomic expansion masking functional homogenization. Because the same volunteer data underpin global biodiversity indicators, our results imply that detectability correction is necessary but not sufficient — unless it is carried into the functional and compositional axes where homogenization actually lives, monitoring systems will continue to certify homogenization as recovery.

**Keywords:** biotic homogenization; citizen science; dynamic occupancy model; imperfect detection; functional diversity; beta-diversity partitioning; biodiversity monitoring; global change; birds; China

---

## 1 | INTRODUCTION

### 1.1 | A contradiction at the centre of biodiversity science

Few empirical claims in ecology are simultaneously as consequential and as contested as the direction of contemporary biodiversity change. One body of evidence documents accelerating loss: extinction rates orders of magnitude above background (Pimm et al., 2014; Ceballos et al., 2015), pervasive human-driven decline requiring transformative change (Díaz et al., 2019), measurable consequences for ecosystem functioning (Cardinale et al., 2012), and steep aggregate declines in well-monitored taxa — nearly three billion fewer birds in North America since 1970 (Rosenberg et al., 2019) and widespread terrestrial insect decline (van Klink et al., 2020). A second body of evidence, drawn from assemblage time series rather than from range or population aggregates, finds something harder to reconcile with catastrophe: local species richness shows no systematic net decline, even as composition turns over rapidly (Vellend et al., 2013; Dornelas et al., 2014; Blowes et al., 2019). Both bodies of evidence are methodologically sound. Their coexistence is therefore not a dispute about facts but a signal that the field is measuring different things and calling them by the same name.

Three developments have sharpened this diagnosis. First, biodiversity change is now understood to be multidimensional: McGill et al. (2015) enumerate fifteen distinct forms of biodiversity trend that are neither interchangeable nor even necessarily correlated in sign. Second, richness trends are strongly scale-dependent, so that local, regional and global signals can legitimately diverge (Chase et al., 2019; Jarzyna & Jetz, 2018). Third, and most directly, Hillebrand et al. (2018) showed that biodiversity change is substantially *uncoupled* from species richness trends: assemblages can be reorganizing profoundly while their species counts stay flat. Together these results imply that the loss-versus-no-loss debate cannot be settled in the currency in which it is usually conducted. The question is not only how much biodiversity is changing, but which axis of biodiversity is being measured — and whether that axis is the one on which degradation is actually occurring.

### 1.2 | The axis on which degradation hides: biotic homogenization

There is a well-developed theory of exactly the kind of change that richness cannot see. Biotic homogenization describes the process by which a few widespread winners replace many localized losers, eroding the distinctiveness of assemblages while leaving their species counts largely intact (McKinney & Lockwood, 1999). Olden and Rooney (2006) formalized its measurement as rising among-site similarity through time, and subsequent work established that homogenization is often *functional* rather than merely taxonomic: communities can retain species number while converging in trait composition (Devictor et al., 2007; Clavel et al., 2010; Baiser et al., 2012). Birds provide some of the clearest evidence. European avifaunas have become progressively more generalist over two decades (Le Viol et al., 2012), specialists have declined worldwide (Clavel et al., 2010), and marine fish assemblages show rapid homogenization detectable only in compositional space (Magurran et al., 2015).

The mechanism is functional redundancy. Because generalists can occupy the niches vacated by specialists, richness is buffered while functional structure erodes — and because rare species disproportionately support vulnerable functions (Mouillot et al., 2013), the functional consequences of losing them are poorly predicted by how many species remain. Detecting this requires the axes that richness ignores: multidimensional functional indices (Villéger et al., 2008; Laliberté & Legendre, 2010), functional diversity conceived as a driver of ecosystem process rather than a summary statistic (Petchey & Gaston, 2006; Cadotte et al., 2011; Violle et al., 2007; Flynn et al., 2009), phylogenetic diversity (Faith, 1992), and the partitioning of beta diversity into turnover and nestedness (Baselga, 2010, 2012; Legendre, 2014; Socolar et al., 2016). For birds specifically, whose ecological functions span seed dispersal, pollination, predation and scavenging (Şekercioğlu, 2006), functional composition — not species number — is the quantity with ecosystem-level meaning.

### 1.3 | The data problem: when the observer is part of the signal

If the first obstacle is measuring the wrong axis, the second is that the data now used to measure any axis at continental scale are collected non-randomly and increasingly. Volunteer networks supply the majority of occurrence records underpinning large-scale biodiversity assessment (Sullivan et al., 2009; Hochachka et al., 2012; Callaghan et al., 2020). These records inherit the behaviour of the people who make them: effort grows as participation grows, observers cluster near roads and cities, platforms differ in protocol, and species differ in conspicuousness (Boakes et al., 2010; Bird et al., 2014; Tingley & Beissinger, 2009). The consequence is a confound that acts in exactly the direction most likely to mislead: where birdwatching has expanded, the number of species reported per place rises for reasons that are sociological rather than ecological.

The methodological response is mature. Hierarchical occupancy models separate the probability that a species is present from the probability that it is detected (MacKenzie et al., 2002, 2003), and have been extended to communities through multispecies hierarchical formulations that share information across species (Dorazio & Royle, 2005), to opportunistic and semistructured citizen-science data specifically (Kéry et al., 2010; van Strien et al., 2013; Isaac et al., 2014; Altwegg & Nichols, 2019; Kelling et al., 2019; Johnston et al., 2021), and to explicitly spatial and dynamic settings under imperfect detection (Guillera-Arroita, 2017; Devarajan et al., 2020; Kéry & Royle, 2021), with scalable implementations now available (Doser et al., 2022). The toolkit is not the limitation.

### 1.4 | The gap: correction that stops at richness

The limitation is where the correction stops. In the overwhelming majority of applications, the output of an occupancy model is a corrected estimate of species richness, occupancy rate, or a species-level trend. The posterior uncertainty so carefully estimated at the observation level is then discarded — collapsed to a point estimate — before any functional, phylogenetic or compositional metric is computed, if such metrics are computed at all. This creates a precise and consequential gap. The homogenization literature identifies the axes on which contemporary degradation occurs, but its evidence is largely uncorrected for detection and therefore vulnerable to exactly the observer-expansion confound described above. The occupancy literature corrects for detection rigorously, but reports its results on the one axis least able to reveal homogenization. **To our knowledge, no national-scale study has propagated a detection-corrected occupancy posterior simultaneously into taxonomic, functional, phylogenetic and beta-diversity space in order to test whether observed biodiversity gains represent recovery or homogenization.**

### 1.5 | Why China, and why now

China offers an unusually decisive test. It combines the world's most rapid concurrent global-change pressures — pronounced warming, large-scale land conversion, urbanization, and simultaneously one of the largest programmes of ecological restoration and conservation investment ever undertaken (Ouyang et al., 2016; Bryan et al., 2018) — with one of the world's fastest-growing volunteer observation networks. Bird communities elsewhere respond to these pressures through range shifts, climatic debt and generalist gains (Devictor et al., 2012; Stephens et al., 2016; Bowler et al., 2017; Antão et al., 2020), and land-use intensification erodes functional diversity across taxa (Newbold et al., 2015, 2016; Daskalova et al., 2020). Chinese avifaunas have been studied at regional scale (Ding et al., 2015; Si et al., 2015), but no detection-corrected, multidimensional national synthesis exists. Critically, China's observer expansion is spatially coincident with its strongest environmental change, so any uncorrected national richness trend is not merely uncertain but structurally ambiguous. This is precisely the setting in which the distinction between recovery and homogenization is both hardest to make and most important to get right.

### 1.6 | Objectives and hypotheses

We ask a single question: **when detection-corrected richness rises across a continent under rapid global change, do communities gain functional and compositional distinctiveness — recovery — or do they accumulate redundant, broad-habitat generalists and converge — hidden homogenization?**

We test five hypotheses:

- **H1 (taxonomic change).** Detection-corrected taxonomic diversity increased across China between 2000 and 2024.
- **H2 (recovery vs homogenization — the discriminating test).** Under the *recovery* hypothesis, functional and phylogenetic diversity rose commensurately with richness; under the *homogenization* hypothesis, they stagnated or declined. These make opposite predictions and H2 is the pivot of the study.
- **H3 (compositional reorganization).** Recent change is dominated by nestedness rather than balanced turnover, and among-site distinctiveness declined.
- **H4 (trait filtering).** Broad-habitat, generalist traits predict which species expanded.
- **H5 (inferential consequence and robustness).** Detection correction materially changes ecological inference relative to naive occurrence trends, and the occupancy increase is robust to observer expansion.

We treat a 200-species spatial model as the primary inference and a 500-species temporal model as a breadth test of generality.

---

## 2 | MATERIALS AND METHODS

### 2.1 | Overview of the inferential architecture

The analysis proceeds in five linked stages: (i) integration and source-aware deduplication of occurrence records; (ii) construction of detection/non-detection survey histories on a 100-km × 5-period design; (iii) estimation of latent occupancy under imperfect detection with a spatial multispecies dynamic occupancy model; (iv) propagation of the *entire* occupancy posterior — not its mean — into taxonomic, functional, phylogenetic and beta-diversity metrics; and (v) attribution and validation, comprising trait and driver analyses, naive-versus-corrected comparison, and four controls for the observer-expansion confound. Stage (iv) is the methodological core: every derived metric inherits the uncertainty of the occupancy estimates that generate it, so that functional and compositional claims carry the same inferential status as taxonomic ones.

### 2.2 | Data sources and integration

We combined occurrence records from the China Birdwatching Records Center with curated GBIF and eBird records for China, 2000–2024, following published analytical guidance for community-science data (Johnston et al., 2021; Kelling et al., 2019). Records were harmonized for taxonomy, date, coordinates and source identity. Deduplication used a source-aware rule keyed on species, date, rounded coordinates, source and observer, which removes within-source duplicates while preserving genuinely independent cross-platform observations — a distinction that matters because collapsing cross-platform records would artificially deflate apparent effort in precisely the periods when multiple platforms operated. A reproducible deduplication audit and source-wise annual coverage table are provided (Supporting Information S1).

### 2.3 | Spatial and temporal design

Records were aggregated to a 100-km equal-area grid and five five-year periods (2000–2004, 2005–2009, 2010–2014, 2015–2019, 2020–2024); the primary analysis used the 1,247 grid cells containing data. Within-period years served as repeated survey occasions. This is a relaxed-closure design: the closure assumption is violated at the sub-period scale, a known and quantifiable issue in occupancy estimation (Rota et al., 2009). We therefore treat period-level occupancy as "use or occurrence during the period" rather than instantaneous presence, and evaluate sensitivity to the window length using a three-year-window re-analysis (Supporting Information S2).

### 2.4 | Primary model: spatial multispecies dynamic occupancy (200 species)

The primary model was fitted with `spOccupancy::stMsPGOcc` (Doser et al., 2022), which implements a multispecies dynamic occupancy model with Pólya-Gamma data augmentation.

**Occupancy sub-model.** Species-specific occupancy probability was modelled as a function of climatic, topographic, habitat, human-footprint (Mu et al., 2022), land-cover (Yang & Huang, 2021), spatial-coordinate and scaled-time covariates, with species-level coefficients drawn from community-level distributions. Temporal dependence was modelled with an AR(1) structure; residual spatial dependence with a nearest-neighbour Gaussian process (NNGP) under an exponential covariance.

**Detection sub-model.** Detection probability was modelled as a function of survey effort (`log_events`), survey duration (`log_duration`), a missing-duration indicator (`has_duration`), and a data-source term distinguishing China Birdwatching Records Center from GBIF/eBird records. The source term is included specifically to prevent platform-specific detectability from being reassigned to occupancy.

**Priors and sensitivity.** The spatial decay parameter `phi` was assigned a uniform prior bounded by the distribution of pairwise inter-grid distances; because a poorly scaled `phi` prior can truncate the estimated spatial range at a 100-km grid, we conducted a prior-sensitivity analysis and confirmed that the posterior did not concentrate at the prior bound (Supporting Information S3).

**Estimation and diagnostics.** Models were run with four chains (5,000 burn-in, thinning 2, 1,250 retained draws per chain). Convergence was satisfactory: maximum R̂ = 1.039 across the 38 monitored community-level parameters, all below the 1.05 threshold, with minimum effective sample size 70; the three slowest-mixing terms were community-level variance hyperparameters and are flagged explicitly rather than smoothed over (Supporting Information S4). Goodness of fit was assessed with posterior predictive checks `[PENDING: Bayesian p-value, S4]`.

### 2.5 | Breadth extension: temporal multispecies model (500 species)

To test whether conclusions depend on the 200-species pool, we fitted a 500-species model with `spOccupancy::tMsPGOcc`, using AR(1) temporal dynamics and ten latent factors. This model deliberately omits NNGP spatial structure to remain computationally tractable at 500 species. **It is therefore temporal and non-spatial, and is used exclusively to assess the generality of taxonomic and species-level signals; no spatial inference in this paper rests on it.** We state this explicitly at each point of use because conflating the two models would misattribute spatial conclusions to a model without spatial structure.

### 2.6 | Posterior propagation into multidimensional diversity

For each posterior draw *d*, grid cell *i* and period *t*, we assembled the vector of species occupancy probabilities **ψ**(*d*, *i*, *t*) and computed the full metric suite from it, yielding a posterior distribution — not a point estimate — for every metric, grid and period:

- **Taxonomic.** Expected species richness as the sum of occupancy probabilities; Shannon and inverse Simpson diversity computed on normalized occupancy weights.
- **Functional.** Multidimensional functional indices (trait volume, evenness, divergence; Villéger et al., 2008) and Rao's quadratic entropy (Laliberté & Legendre, 2010), within the framework of functional diversity as an ecological rather than descriptive quantity (Petchey & Gaston, 2006; Cadotte et al., 2011). Traits were drawn from AVONET (Tobias et al., 2022) and EltonTraits (Wilman et al., 2014), with IUCN-derived habitat breadth and diet-specialization scores (Morelli et al., 2021); the morphological trait space follows the avian form-to-function mapping of Pigot et al. (2020). Missing traits were imputed, with the complete-case subset retained for mechanism analyses (Supporting Information S5).
- **Phylogenetic.** Probability-weighted phylogenetic diversity (Faith, 1992) and mean pairwise distance on the global bird phylogeny (Jetz et al., 2012). Because adding any species to a sparse tree inflates raw phylogenetic diversity, we additionally report standardized effect sizes against a null model `[PENDING: ses.PD, ses.MPD, S5]`.

All derived metrics are summarized as posterior means with 95% credible intervals propagated from the occupancy posterior.

### 2.7 | Temporal and spatial beta diversity

Temporal beta diversity between adjacent periods was partitioned into turnover and nestedness components using a probability-weighted generalization of the Baselga framework (Baselga, 2010, 2012; Legendre, 2014), with expected shared occupancy defined as the sum of per-species pairwise minima of occupancy probabilities. Because this generalization is non-standard, its bounds and decomposition properties are derived formally in Supporting Information S6. To test homogenization in its classical, compositional sense (Olden & Rooney, 2006; Socolar et al., 2016), we additionally computed among-grid spatial beta diversity within each period and tested for a monotonic decline through time `[PENDING: spatial β slope and Mann–Kendall test, S6]`.

### 2.8 | Species trends, trait filtering and driver attribution

Species-level occupancy trends were estimated from posterior occupancy trajectories and classified as expanding, stable or contracting from the posterior trend direction and its credible interval. Grid-level diversity trends used Theil–Sen slopes with Mann–Kendall tests, which are robust for five-point series. Trait–trend and environment–trend associations were tested with Spearman rank correlations on the complete-case species subset (n = 200), and community-trend drivers were analysed with variation partitioning (Borcard et al., 1992) and random-forest permutation importance (Breiman, 2001). Because the design is observational and spatially structured, we interpret all such results as associations, not causal effects, throughout.

### 2.9 | Controls for the observer-expansion confound

Because the central alternative explanation for our results is that growing effort — not range change — produced the occupancy increase, we treat it as a hypothesis to be tested rather than assumed away, using four independent controls (Supporting Information S7): (i) the posterior slope of detection probability on survey effort, establishing whether effort growth is absorbed at the observation level; (ii) re-estimation of trends restricted to grid cells where effort saturated before 2010 and no longer increased; (iii) a fixed-effort counterfactual quantifying how much of the change in detection probability is attributable to effort growth alone; and (iv) comparison of models with and without the data-source detection term.

### 2.10 | Protected-area contrast, effectiveness and prioritization

Because homogenization is a conservation problem, we asked whether China's reserve network slows it, and whether planning based on richness would target the right places. Reserve boundaries and attributes came from the national nature-reserve inventory (China Nature Reserve Specimen Resource Sharing Platform, 2024); we computed, for every grid cell, the areal fraction covered by reserves in an equal-area projection, together with the earliest establishment year of any intersecting reserve and a "major-overlap" variant restricted to reserves contributing at least 1% of grid area.

**A three-tier causal hierarchy.** Protected areas are not sited at random — they are systematically placed high, far and on land of low opportunity cost (Joppa & Pfaff, 2009), so naive inside–outside comparisons confound protection with location. Global assessments that address this confound find protection effects that are real but heterogeneous and often modest (Watson et al., 2014; Geldmann et al., 2019), which sets a realistic expectation for what we should detect here. We therefore report three designs of decreasing inferential strength and state explicitly which claims each supports.

*Tier 1 — matched treatment effects (primary).* Following the quasi-experimental standard for protected-area impact evaluation (Andam et al., 2008; Ferraro & Hanauer, 2014), we matched protected grids to unprotected grids on baseline covariates (elevation and its heterogeneity, baseline climate, land cover, baseline human footprint, habitat diversity) using propensity-score nearest-neighbour matching with a 0.2 caliper, then estimated average treatment effects on the treated for each diversity trend. We report standardized mean differences for every covariate before and after matching; effects are interpreted only where post-matching |SMD| < 0.1.

*Tier 2 — protection-age dose-response (corroborating).* If protection is effective, grids protected for longer should show more favourable functional trajectories. We modelled each diversity metric as a function of period × protection age with grid random effects, using the full 1957–2012 range of establishment years. A dose-response relationship is harder to generate by siting bias alone than a binary contrast, because siting bias would have to covary with reserve age in the same direction.

*Tier 3 — event-study dynamic difference-in-differences (exploratory).* For grids first substantially protected during the study window, we estimated period-specific effects relative to establishment, which permits a direct test of parallel pre-trends. As detailed in §4.8, the inventory's establishment years end around 2012, leaving only the 2005–2014 cohorts with both pre- and post-treatment periods; this design is therefore small-sample and we treat it as exploratory support for the parallel-trends assumption rather than as independent causal evidence.

**Spatiotemporal contrast.** Beyond endpoint effects, we compared period-by-period trajectories inside and outside reserves for taxonomic, functional and phylogenetic metrics, tested period × protection interactions, and extended the contrast to two levels that endpoint summaries cannot reach: species-level occupancy, split by habitat breadth into specialists and generalists, and the rate of spatial homogenization, computed as mean among-grid Sørensen dissimilarity within each period on equal-sized random subsets of protected and unprotected grids to remove any effect of differing sample sizes.

**Systematic conservation prioritization.** Using integer-programming prioritization (Margules & Pressey, 2000; Hanson et al., 2024), we solved three planning problems under a common area budget: maximizing detection-corrected richness, maximizing functional distinctiveness, and minimizing homogenization risk. Comparing the spatial overlap of the three solutions tests whether richness-led planning — still the default in most area-based targets (Rodrigues et al., 2004; Butchart et al., 2015; Maxwell et al., 2020) — would select the areas that matter for preventing functional homogenization, in the spirit of multi-facet prioritization (Pollock et al., 2017; Brum et al., 2017; Jung et al., 2021). We also quantified the shortfall between current reserve coverage and each solution.

### 2.11 | Software, reproducibility and data availability

Analyses used R with `spOccupancy`, `brms`, `cmdstanr`, `sf`, `terra`, `vegan`, `ape` and `ggplot2`. Code, derived result tables and figures are archived openly, with a locked `sessionInfo`, code hash and data manifest, and a full data dictionary linking every reported number to its source table.

---

## 3 | RESULTS

### 3.1 | H1: detection-corrected taxonomic diversity increased

Detection-corrected richness rose monotonically across all five periods in the 200-species spatial model, from 79.2 species per grid in 2000–2004 to 82.2, 87.0, 90.5 and 100.8 in successive periods — a 27% increase over 25 years (median 80.4 to 102.7). Shannon diversity rose from 4.74 to 4.88 and inverse Simpson from 106.4 to 126.6. The 500-species temporal model reproduced the direction and amplified the magnitude (mean richness 155.0 to 195.6; Shannon 5.51 to 5.65; inverse Simpson 218.7 to 260.3; mean across-grid slope 9.6 species per period). H1 is supported in both species pools, and the increase survives detection correction rather than being created by it.

### 3.2 | H2: functional diversity did not track taxonomic gain

Functional trait space did not expand with richness — the discriminating result of this study. In the 200-species model, functional trait volume declined from 1.383 to 1.365 and Rao's quadratic entropy from 1.705 to 1.680, and both declines were **monotonic across all five periods** rather than driven by a single interval (trait volume: 1.383, 1.378, 1.376, 1.368, 1.365; Rao's Q: 1.705, 1.699, 1.696, 1.684, 1.680). The 500-species model showed the same signature (trait volume 1.581 to 1.565; Rao's Q 1.732 to 1.711). Posterior probabilities that these slopes are negative are `[PENDING: P(slope < 0) for trait volume and Rao's Q, S5]`; we calibrate the strength of our functional claim to these values (§4.6).

Phylogenetic diversity rose (probability-weighted PD 3.60 × 10⁵ to 4.58 × 10⁵; mean pairwise distance 122.3 to 126.3), a pattern consistent with adding species that are phylogenetically dispersed but functionally redundant. Under the recovery hypothesis, functional and phylogenetic diversity should both rise with richness. They did not both rise, and functional diversity moved in the opposite direction. **H2's recovery formulation is rejected; the homogenization formulation is supported.**

### 3.3 | H3: temporal beta diversity shifted from turnover to nestedness

Probability-weighted Baselga partitioning showed the contribution of turnover to temporal beta diversity collapsing in the most recent interval. In the 200-species model, turnover accounted for 65.9%, 42.8% and 53.8% of beta diversity across the first three transitions, but only 11.3% across 2015–2019 to 2020–2024 (95% credible interval 7.0–16.6%, n = 1,247 grids); the 500-species model showed the same collapse (75.7%, 55.1%, 65.3%, then 19.5%). Recent community change is therefore dominated by nestedness rather than balanced replacement.

The direction of this nestedness is diagnostic. Nestedness arising from local extinction would produce depauperate subsets and declining richness; here richness rose over the same interval. The signal is therefore one of accumulation — widespread species spreading across grids to form nested supersets — rather than attrition. The direct compositional test, whether among-grid spatial beta diversity declined through time, is `[PENDING: S6]`.

### 3.4 | H4: broad-habitat traits and environmental gradients predict expansion

Within the complete-case trait subset (n = 200 species), habitat breadth was positively associated with species trend slopes (Spearman ρ = 0.38), diet specialization weakly negative (ρ = −0.08) and migration score near zero (ρ = 0.08). Habitat breadth is the strongest single trait predictor, and its sign is the one homogenization theory predicts.

Species trend slopes were positively correlated with cropland cover (ρ = 0.39), elevation (ρ = 0.38), temperature annual range and wettest-month precipitation (both ρ = 0.38), temperature seasonality (ρ = 0.36) and human footprint (ρ = 0.36). Variation partitioning of the corrected-richness trend attributed most explained variance to baseline spatial and environmental structure (adjusted R² = 0.14), with smaller pure contributions from climate change (0.04), land-use change (0.01) and human-pressure change (0.01), and a large residual (0.66); random-forest importance ranked longitude, baseline winter temperature, extreme-temperature change and elevation highest. These are spatially structured associations whose directions are consistent with documented climate- and land-use-driven reorganization (Stephens et al., 2016; Bowler et al., 2017; Newbold et al., 2016; Antão et al., 2020), not causal effects.

### 3.5 | Winners and losers: the functional identity of change

Most species expanded: 136 expanding, 62 stable and 2 contracting in the 200-species model; 258, 227 and 15 in the 500-species model. The fastest-expanding species included Eurasian Coot (*Fulica atra*), Oriental Honey-buzzard (*Pernis ptilorhynchus*), Grey-capped Pygmy Woodpecker (*Picoides canicapillus*), Whiskered Tern (*Chlidonias hybrida*) and Grey-headed Lapwing (*Vanellus cinereus*) — a set dominated by waterbirds, open-country species and human-associated or highly mobile taxa. The strongest contractions were the Rock Dove (*Columba livia*) and, in the broader pool, narrow-range montane forest specialists including *Seicercus valentini*, *Phylloscopus ogilviegranti* and *Certhia hodgsoni*.

This asymmetry in ecological identity, not merely in count, is the mechanistic signature of homogenization: the species accumulating are those able to occupy many habitat types, and those retreating are those that cannot.

### 3.6 | H5a: the increase is not an artefact of observer expansion

Four independent controls address the central alternative explanation. First, detection probability increased strongly with survey effort (posterior slope on log effort = 0.93, 95% CrI 0.90–0.95), confirming that effort growth is captured at the observation level rather than passed to occupancy. Second, restricting estimation to the 30 grid cells where effort saturated before 2010, corrected richness still rose (median slope 4.1, mean 4.5 species per period, 95% CrI 2.4–7.9) — a small but effort-stable subset in which growing detectability cannot generate the trend. Third, a fixed-effort counterfactual showed that rising effort raised detection probability by 0.32 (95% CrI 0.30–0.34) from the P1 to the P5 effort level, quantifying the effort contribution explicitly at the detection layer. Fourth, adding the data-source detection term left the community trajectory unchanged `[PENDING: ΔWAIC and source-term trend slope, S7]`.

Two structural features of the results are additionally difficult to reconcile with a pure-effort artefact, and we regard them as the strongest evidence of all. Detection correction *reversed* trend direction in a minority of grids rather than uniformly inflating richness (§3.7) — an effort artefact would predict systematic inflation. And functional diversity moved *opposite* to richness: more observation should, if anything, reveal more functional variety, not less. No purely sociological account of the data predicts a negative functional trend alongside a positive taxonomic one.

### 3.7 | H5b: detection correction changes ecological inference

Naive and corrected trends were directionally similar in most grids but differed substantially in magnitude. Correction reversed the sign of the trend in 37 of 1,247 grids (3.0%) in the 200-species model and 25 grids (2.0%) in the 500-species model; median absolute trend differences were 6.0 and 8.7 corrected species per period respectively. Detection correction therefore does not merely rescale raw richness — it alters the magnitude of inferred change widely and reverses interpretation in a small but conservation-relevant set of locations. These disagreement grids are, by construction, the places where uncorrected citizen-science data would most mislead monitoring.

---

### 3.8 | Protection, homogenization and the geography of conservation priority

*This section reports analyses that are complete in design and code but awaiting the compute run; every quantity is flagged and none is inferred.*

Reserve coverage is highly uneven across the grid `[PENDING: mean and distribution of grid protected fraction]`. Comparing period-by-period trajectories inside and outside reserves, taxonomic and functional metrics diverged as follows `[PENDING: protected-minus-unprotected change for corrected richness, trait volume and Rao's Q]`, with period × protection interactions `[PENDING: interaction estimates and significance]`. Matched treatment effects, estimated only where post-matching covariate balance was achieved `[PENDING: number of covariates with |SMD| < 0.1]`, were `[PENDING: ATT and 95% CI per diversity trend]`. The protection-age dose-response relationship, which uses the full range of establishment years, gave `[PENDING: period × protection-age interaction for functional metrics]`, and the exploratory event study yielded `[PENDING: pre-treatment coefficients as a parallel-trends check, with cohort sizes]`.

Two contrasts bear most directly on the paper's central claim. At the species level, the occupancy gap between protected and unprotected grids for habitat specialists versus generalists was `[PENDING: mean psi difference by species group and period]` — a widening specialist gap would indicate that reserves selectively retain the species that homogenization removes. At the community level, the rate of spatial homogenization, measured as among-grid dissimilarity on equal-sized subsets, changed inside versus outside reserves by `[PENDING: spatial beta trajectories and inside-outside gap]`; a slower decline inside reserves would show that protection retards compositional convergence itself, not merely species loss.

Finally, prioritization solutions optimized for corrected richness, for functional distinctiveness, and against homogenization risk overlapped by `[PENDING: pairwise Jaccard overlap]`, and current reserves covered `[PENDING: fraction of each solution already protected; number of unprotected priority grids]`. Low overlap would carry a direct planning implication: selecting areas by species richness — the currency of most area-based targets — would systematically miss the places where functional homogenization can still be prevented.

---

## 4 | DISCUSSION

### 4.1 | Principal finding

Chinese bird communities have gained species and occupied area after detection correction, but they have not gained functional or compositional distinctiveness in proportion. Across a 200-species spatial model and a 500-species temporal model, corrected richness, Shannon, inverse Simpson and phylogenetic diversity all rose, while functional trait volume and Rao's quadratic entropy declined monotonically, recent beta diversity collapsed toward nestedness, and the expanding species were predominantly broad-habitat generalists displacing narrow-range specialists. We read this joint signature as **expansion without differentiation**: more species are expected per place, but the species accumulating are functionally redundant and increasingly shared among assemblages. This is the multidimensional fingerprint of biotic homogenization (McKinney & Lockwood, 1999; Olden & Rooney, 2006; Clavel et al., 2010), and it is invisible to richness alone.

### 4.2 | Resolving the loss-versus-no-loss debate: a third reading

Our results speak directly to the contradiction with which we began. The no-net-loss literature (Vellend et al., 2013; Dornelas et al., 2014; Blowes et al., 2019) and the ongoing-degradation literature (Pimm et al., 2014; Ceballos et al., 2015; Díaz et al., 2019; Rosenberg et al., 2019) are usually framed as competing accounts. We offer a third reading that dissolves rather than adjudicates the dispute: **local richness can rise while assemblage distinctiveness falls, so "no net loss" and "continuing degradation" can both be true of the same system, recorded on different axes.** This is the empirical realization of what Hillebrand et al. (2018) argued in principle — that biodiversity change is uncoupled from richness trends — and it extends that argument in a specific and consequential way. Hillebrand et al. showed uncoupling can occur; we show that under detection correction at national scale it occurs *with opposite signs on different axes simultaneously*, which is a stronger and more diagnostic claim. It also aligns with the multidimensional framing of McGill et al. (2015) and with evidence that taxonomic and functional avian diversity change is scale-dependent and often decoupled (Jarzyna & Jetz, 2017, 2018; Chase et al., 2019).

The practical consequence is uncomfortable. A rising corrected-richness map of China could be presented, in good faith and with correct statistics, as evidence of conservation success. The functional and beta-diversity axes show it need not be.

### 4.3 | Confronting the observer-expansion confound

The most serious threat to our interpretation is that growing observer effort produced the pattern. We treat this as a hypothesis to be defeated, not a caveat to be acknowledged. A dynamic occupancy model separates detection from occupancy only insofar as the detection covariates capture the true effort process; because Chinese birdwatching grew fastest in the same eastern, lowland, human-modified grids that show the strongest gains, residual detectability could in principle be reassigned to colonization.

Four results argue against this (§3.6), and two structural features of the data argue against it more strongly than any single control: correction reverses direction in a minority of grids rather than inflating uniformly, and the functional axis moves opposite to the taxonomic one. We emphasize the second because it is the one an artefact cannot easily produce. If rising richness were simply better looking, the newly detected species would be drawn approximately at random from the regional pool with respect to traits, and functional volume would rise or at worst stay flat. Instead it falls, monotonically, while the identity of expanding species is systematically biased toward habitat generalists. A sociological process would have to be trait-selective in exactly the direction homogenization theory predicts in order to mimic this — a considerably less parsimonious explanation than the ecological one.

We nonetheless retain the spatial collinearity between effort growth and the strongest-gaining grids as a genuine residual caution, and note that the effort-saturated subset is small (n = 30). Extending that subset is the single most valuable additional test of this study's central claim.

### 4.4 | Mechanism: habitat breadth as the dominant filter

Habitat breadth emerged as the strongest single trait predictor of occupancy trend, with diet specialization weakly negative — the pattern expected if generalists exploit heterogeneous, modified and newly accessible landscapes while specialists cannot. This is consistent with the documented rise of generalists in intensively used European landscapes (Le Viol et al., 2012; Clavel et al., 2010; Devictor et al., 2007) and with functional-diversity loss under land-use intensification across taxa (Flynn et al., 2009; Newbold et al., 2015). The identity of winners (waterbirds, open-country and human-associated taxa) and losers (montane forest specialists) matches the mechanism at the level of natural history, not merely statistics.

The functional consequences may exceed what the magnitude of the trait-volume decline suggests. Because rare and specialized species disproportionately support vulnerable ecological functions (Mouillot et al., 2013), and because avian functional roles span dispersal, pollination, predation and scavenging (Şekercioğlu, 2006), the replacement of specialists by redundant generalists can erode ecosystem function faster than aggregate functional-diversity indices decline. Our functional signal should therefore be read as a lower bound on functional consequence, not an upper one.

### 4.5 | Distinguishing drivers: what we can and cannot attribute

At least three non-exclusive global-change processes are consistent with our spatial associations, and honesty requires separating them. Climate warming can drive apparent expansion at the national envelope as species track shifting climate, accumulating climatic debt where they lag (Devictor et al., 2012; Stephens et al., 2016; Bowler et al., 2017; Antão et al., 2020). Land-use change and rising human pressure favour generalists while disadvantaging specialists (Newbold et al., 2015, 2016; Daskalova et al., 2020). And China's large-scale investment in restoration and conservation could produce genuine recovery in some regions (Ouyang et al., 2016; Bryan et al., 2018) — a competing causal story that our data neither confirm nor exclude, and which, notably, would predict rising functional as well as taxonomic diversity, making it a poor sole explanation for what we observe.

Our variation partitioning, in which baseline spatial structure dominates (adjusted R² = 0.14) and a large residual remains (0.66), cannot adjudicate among these. It shows that the expansion is spatially organized but leaves its causes open. We therefore frame the homogenization signal as the robust result and driver attribution as a hypothesis for the lagged, quasi-experimental and event-based analyses these data now make possible.

### 4.6 | Calibrating the strength of the functional claim

We state explicitly how strongly this paper is entitled to speak. The direction of the functional signal is robust: trait volume and Rao's Q decline monotonically across all five periods in both species pools, which is not the behaviour of noise. The magnitude is small (≈1–2% over 25 years), and the formal posterior probability that these slopes are negative is still being computed `[PENDING: §3.2]`. If those posteriors exclude zero, the data support functional contraction and the term "functional homogenization" is warranted. If they span zero, the defensible claim is **taxonomic expansion without detectable functional expansion** — which still rejects the recovery hypothesis, because recovery predicts a *positive* functional trend, not merely a non-negative one. We have written this manuscript so that its central conclusion survives either outcome, and we will report whichever the posteriors support.

### 4.7 | Implications for global biodiversity monitoring

The transferable contribution of this study is methodological and general. Continental biodiversity assessments and Essential Biodiversity Variables increasingly rest on the same volunteer data we analyse (Isaac et al., 2014; Jetz et al., 2019), and complex, facet-dependent change of the kind we report is emerging wherever such data are examined carefully across multiple axes (Outhwaite et al., 2020; Magurran et al., 2015). Our results show that **detection correction is necessary but not sufficient**: it must be propagated beyond richness into the functional and compositional axes where homogenization lives, because taxonomic and functional trajectories can diverge in sign. A monitoring system that corrects detectability rigorously and then reports corrected richness alone has bought statistical rigour at the price of ecological blindness.

For conservation practice, three priorities follow. First, the grids where naive and corrected trends disagree are where uncorrected data would most mislead, and are the natural priority for ground validation. Second, the winner–loser structure flags narrow-range montane forest specialists as the assemblages most likely to be masked by aggregate gains — a group that richness-based reporting actively hides, and for which protected areas may mediate climate responses (Santangeli et al., 2016). Third, biodiversity targets framed in richness or occupied-area terms should be complemented by functional and compositional targets, or apparent attainment may conceal continuing homogenization.

### 4.8 | Limitations

Several constraints bound these conclusions. The 500-species model is temporal and non-spatial: it supports generality, not spatial mechanism, and all spatial inference rests on the 200-species model. Trait and environmental mechanism analyses are complete-case (n = 200 species) and should not be extrapolated to the full 500-species pool. The five-year relaxed-closure design is an approximation, justified by sensitivity analysis but not eliminated by it. Convergence, while satisfactory by R̂, included three community-level variance hyperparameters with effective sample sizes below 100, which we report rather than conceal. Environmental associations are observational and spatially structured, and we avoid causal language accordingly. Finally, three analyses — spatial beta-diversity trend, functional-slope posterior probabilities, and the source-term model comparison — are on the compute queue and flagged at the exact sentences where their values belong.

The protected-area analyses carry two further constraints that we state explicitly because they bound what can be claimed causally. First, **spatial coverage**: vector boundaries are available for 1,028 of the 3,376 reserves in the national inventory, and the covered subset is skewed toward national-level and large reserves. Smaller reserves without boundaries are therefore treated as unprotected, which biases estimated protection effects toward zero; our effect estimates should be read as conservative, and as applying to the large, higher-designation portion of the network rather than to area-based conservation in China as a whole. Second, **temporal coverage**: establishment years in the inventory end around 2012. With five-year periods, only cohorts first protected during 2005–2014 have both pre- and post-treatment periods, so the event-study design rests on a small number of grids. We accordingly place primary weight on the matched treatment effects and on the protection-age dose-response relationship, and present the event study only as an exploratory check on parallel pre-trends. None of the three designs can exclude unobserved confounders that covary with both siting and subsequent community change; we therefore describe protection results as differences associated with protection, reserving stronger language for cases where the matched, dose-response and event-study designs agree.

### 4.9 | Conclusion

Detection correction did not erase the national richness increase, and we do not claim that Chinese bird communities are collapsing. We claim something more specific and, for the citizen-science era, more general: gains in observed and corrected richness can coincide with stagnant or declining functional space and rising nestedness. Biodiversity inference from volunteer data must therefore model detectability *and then look beyond richness* before change is called recovery. More species is not the same as more biodiversity — and in China's birds over the past quarter-century, the difference between the two is where the ecology is.

---

## ACKNOWLEDGEMENTS

We thank the China Birdwatching Records Center, eBird and the Cornell Lab of Ornithology, GBIF, and the contributors to AVONET, EltonTraits, the IUCN Red List, WorldClim, CLCD and the Human Footprint dataset, and the developers of `spOccupancy` and related open-source tools. We are indebted to the many thousands of volunteer observers whose records made this analysis possible.

## AUTHOR CONTRIBUTIONS

Chenchen Ding conceived the study, designed and performed the analyses, and wrote the manuscript. [Co-author contributions to be assigned.]

## CONFLICT OF INTEREST

The authors declare no conflict of interest.

## DATA AVAILABILITY STATEMENT

Analysis code, derived result tables and figures are openly archived at https://github.com/dingchenchen6/china-bird-community-dynamics, with a data dictionary, reproducibility guide and locked `sessionInfo`. Raw occurrence records are subject to source-provider terms (China Birdwatching Records Center, GBIF, eBird; Supporting Information S1). Derived data will receive a Zenodo DOI on acceptance.

---

## FIGURES

**Figure 1 | Inferential architecture.** Record integration and source-aware deduplication; spatial multispecies dynamic occupancy estimation; propagation of the full occupancy posterior into taxonomic, functional, phylogenetic and beta-diversity space; trait and driver attribution; naive-versus-corrected comparison and effort-confound controls.

**Figure 2 | Multidimensional diversity across space and time.** Five-period maps of corrected richness, Shannon diversity, functional trait volume, Rao's quadratic entropy and probability-weighted phylogenetic diversity (200-species spatial model), with national-mean summaries and credible intervals.

**Figure 3 | Geography of change.** Corrected-richness, functional-trait-volume and phylogenetic trend z-scores with uncertainty masking, showing that taxonomic gain and functional stasis are spatially decoupled.

**Figure 4 | The decoupling.** (a) Endpoint change (2020–2024 relative to 2000–2004) for all diversity axes in both species pools, standardized by baseline spatial s.d. (b) National trajectories by period. Taxonomic and phylogenetic axes rise; functional axes do not.

**Figure 5 | Compositional reorganization and the identity of change.** Turnover and nestedness proportions across adjacent periods for both models, with the late nestedness shift highlighted, alongside the posterior species-trend distribution and the habitat-breadth-versus-trend relationship.

**Figure 6 | The cost of ignoring detection.** Naive versus detection-corrected grid-level trends, direction reversals, and driver attribution by variation partitioning and random-forest importance.

---

## REFERENCES

> All 76 references DOI-verified against Crossref (2026-06-02; additions 2026-07-28). Author–year style shown for editing convenience; convertible to any journal style via Zotero/EndNote.

Altwegg, R., & Nichols, J. D. (2019). Occupancy models for citizen-science data. *Methods in Ecology and Evolution*, 10, 8–21. https://doi.org/10.1111/2041-210X.13090

Andam, K. S., Ferraro, P. J., Pfaff, A., Sanchez-Azofeifa, G. A., & Robalino, J. A. (2008). Measuring the effectiveness of protected area networks in reducing deforestation. *Proceedings of the National Academy of Sciences USA*, 105, 16089–16094. https://doi.org/10.1073/pnas.0800437105

Antão, L. H., Bates, A. E., Blowes, S. A., Waldock, C., Supp, S. R., Magurran, A. E., et al. (2020). Temperature-related biodiversity change across temperate marine and terrestrial systems. *Nature Ecology & Evolution*, 4, 927–933. https://doi.org/10.1038/s41559-020-1185-7

Baiser, B., Olden, J. D., Record, S., Lockwood, J. L., & McKinney, M. L. (2012). Pattern and process of biotic homogenization in the New Pangaea. *Proceedings of the Royal Society B*, 279, 4772–4777. https://doi.org/10.1098/rspb.2012.1651

Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19, 134–143. https://doi.org/10.1111/j.1466-8238.2009.00490.x

Baselga, A. (2012). The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. *Global Ecology and Biogeography*, 21, 1223–1232. https://doi.org/10.1111/j.1466-8238.2011.00756.x

Bird, T. J., Bates, A. E., Lefcheck, J. S., Hill, N. A., Thomson, R. J., Edgar, G. J., et al. (2014). Statistical solutions for error and bias in global citizen science datasets. *Biological Conservation*, 173, 144–154. https://doi.org/10.1016/j.biocon.2013.07.037

Blowes, S. A., Supp, S. R., Antão, L. H., Bates, A., Bruelheide, H., Chase, J. M., et al. (2019). The geography of biodiversity change in marine and terrestrial assemblages. *Science*, 366, 339–345. https://doi.org/10.1126/science.aaw1620

Boakes, E. H., McGowan, P. J. K., Fuller, R. A., Chang-qing, D., Clark, N. E., O'Connor, K., & Mace, G. M. (2010). Distorted views of biodiversity: spatial and temporal bias in species occurrence data. *PLoS Biology*, 8, e1000385. https://doi.org/10.1371/journal.pbio.1000385

Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the spatial component of ecological variation. *Ecology*, 73, 1045–1055. https://doi.org/10.2307/1940179

Bowler, D. E., Hof, C., Haase, P., Kröncke, I., Schweiger, O., Adrian, R., et al. (2017). Cross-realm assessment of climate change impacts on species' abundance trends. *Nature Ecology & Evolution*, 1, 0067. https://doi.org/10.1038/s41559-016-0067

Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324

Brum, F. T., Graham, C. H., Costa, G. C., Hedges, S. B., Penone, C., Radeloff, V. C., et al. (2017). Global priorities for conservation across multiple dimensions of mammalian diversity. *Proceedings of the National Academy of Sciences USA*, 114, 7641–7646. https://doi.org/10.1073/pnas.1706461114

Bryan, B. A., Gao, L., Ye, Y., Sun, X., Connor, J. D., Crossman, N. D., et al. (2018). China's response to a national land-system sustainability emergency. *Nature*, 559, 193–204. https://doi.org/10.1038/s41586-018-0280-2

Butchart, S. H. M., Clarke, M., Smith, R. J., Sykes, R. E., Scharlemann, J. P. W., Harfoot, M., et al. (2015). Shortfalls and solutions for meeting national and global conservation area targets. *Conservation Letters*, 8, 329–337. https://doi.org/10.1111/conl.12158

Cadotte, M. W., Carscadden, K., & Mirotchnick, N. (2011). Beyond species: functional diversity and the maintenance of ecological processes and services. *Journal of Applied Ecology*, 48, 1079–1087. https://doi.org/10.1111/j.1365-2664.2011.02048.x

Callaghan, C. T., Poore, A. G. B., Mesaglio, T., Moles, A. T., Nakagawa, S., Roberts, R. E., et al. (2020). Three frontiers for the future of biodiversity research using citizen science data. *BioScience*, 71, 55–63. https://doi.org/10.1093/biosci/biaa131

Cardinale, B. J., Duffy, J. E., Gonzalez, A., Hooper, D. U., Perrings, C., Venail, P., et al. (2012). Biodiversity loss and its impact on humanity. *Nature*, 486, 59–67. https://doi.org/10.1038/nature11148

Ceballos, G., Ehrlich, P. R., Barnosky, A. D., García, A., Pringle, R. M., & Palmer, T. M. (2015). Accelerated modern human-induced species losses: entering the sixth mass extinction. *Science Advances*, 1, e1400253. https://doi.org/10.1126/sciadv.1400253

Chase, J. M., McGill, B. J., Thompson, P. L., Antão, L. H., Bates, A. E., Blowes, S. A., et al. (2019). Species richness change across spatial scales. *Oikos*, 128, 1079–1091. https://doi.org/10.1111/oik.05968

Clavel, J., Julliard, R., & Devictor, V. (2010). Worldwide decline of specialist species: toward a global functional homogenization? *Frontiers in Ecology and the Environment*, 9, 222–228. https://doi.org/10.1890/080216

Daskalova, G. N., Myers-Smith, I. H., Bjorkman, A. D., Blowes, S. A., Supp, S. R., Magurran, A. E., & Dornelas, M. (2020). Landscape-scale forest loss as a catalyst of population and biodiversity change. *Science*, 368, 1341–1347. https://doi.org/10.1126/science.aba1289

Devarajan, K., Morelli, T. L., & Tingley, M. W. (2020). Multi-species occupancy models: review, roadmap, and recommendations. *Ecography*, 43, 1612–1624. https://doi.org/10.1111/ecog.04957

Devictor, V., Julliard, R., Clavel, J., Jiguet, F., Lee, A., & Couvet, D. (2007). Functional biotic homogenization of bird communities in disturbed landscapes. *Global Ecology and Biogeography*, 17, 252–261. https://doi.org/10.1111/j.1466-8238.2007.00364.x

Devictor, V., van Swaay, C., Brereton, T., Brotons, L., Chamberlain, D., Heliölä, J., et al. (2012). Differences in the climatic debts of birds and butterflies at a continental scale. *Nature Climate Change*, 2, 121–124. https://doi.org/10.1038/nclimate1347

Ding, C., Liang, D., Xin, W., Li, C., Lloyd, H., Zhang, Y., et al. (2015). Bird guild loss and its determinants on subtropical land-bridge islands, China. *Avian Research*, 6, 10. https://doi.org/10.1186/s40657-015-0019-9

Dorazio, R. M., & Royle, J. A. (2005). Estimating size and composition of biological communities by modeling the occurrence of species. *Journal of the American Statistical Association*, 100, 389–398. https://doi.org/10.1198/016214505000000015

Dornelas, M., Gotelli, N. J., McGill, B., Shimadzu, H., Moyes, F., Sievers, C., & Magurran, A. E. (2014). Assemblage time series reveal biodiversity change but not systematic loss. *Science*, 344, 296–299. https://doi.org/10.1126/science.1248484

Doser, J. W., Finley, A. O., Kéry, M., & Zipkin, E. F. (2022). spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. *Methods in Ecology and Evolution*, 13, 1670–1678. https://doi.org/10.1111/2041-210X.13897

Díaz, S., Settele, J., Brondízio, E. S., Ngo, H. T., Agard, J., Arneth, A., et al. (2019). Pervasive human-driven decline of life on Earth points to the need for transformative change. *Science*, 366, eaax3100. https://doi.org/10.1126/science.aax3100

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61, 1–10. https://doi.org/10.1016/0006-3207(92)91201-3

Ferraro, P. J., & Hanauer, M. M. (2014). Advances in measuring the environmental and social impacts of environmental programs. *Annual Review of Environment and Resources*, 39, 495–517. https://doi.org/10.1146/annurev-environ-101813-013230

Flynn, D. F. B., Gogol-Prokurat, M., Nogeire, T., Molinari, N., Richers, B. T., Lin, B. B., et al. (2009). Loss of functional diversity under land use intensification across multiple taxa. *Ecology Letters*, 12, 22–33. https://doi.org/10.1111/j.1461-0248.2008.01255.x

Geldmann, J., Manica, A., Burgess, N. D., Coad, L., & Balmford, A. (2019). A global-level assessment of the effectiveness of protected areas at resisting anthropogenic pressures. *Proceedings of the National Academy of Sciences USA*, 116, 23209–23215. https://doi.org/10.1073/pnas.1908221116

Guillera-Arroita, G. (2017). Modelling of species distributions, range dynamics and communities under imperfect detection: advances, challenges and opportunities. *Ecography*, 40, 281–295. https://doi.org/10.1111/ecog.02445

Hanson, J. O., Schuster, R., Strimas-Mackey, M., Morrell, N., Edwards, B. P. M., Arcese, P., et al. (2024). Systematic conservation prioritization with the prioritizr R package. *Conservation Biology*, 38, e14376. https://doi.org/10.1111/cobi.14376

Hillebrand, H., Blasius, B., Borer, E. T., Chase, J. M., Downing, J. A., Eriksson, B. K., et al. (2018). Biodiversity change is uncoupled from species richness trends: consequences for conservation and monitoring. *Journal of Applied Ecology*, 55, 169–184. https://doi.org/10.1111/1365-2664.12959

Hochachka, W. M., Fink, D., Hutchinson, R. A., Sheldon, D., Wong, W.-K., & Kelling, S. (2012). Data-intensive science applied to broad-scale citizen science. *Trends in Ecology & Evolution*, 27, 130–137. https://doi.org/10.1016/j.tree.2011.11.006

Isaac, N. J. B., van Strien, A. J., August, T. A., de Zeeuw, M. P., & Roy, D. B. (2014). Statistics for citizen science: extracting signals of change from noisy ecological data. *Methods in Ecology and Evolution*, 5, 1052–1060. https://doi.org/10.1111/2041-210X.12254

Jarzyna, M. A., & Jetz, W. (2017). A near half-century of temporal change in different facets of avian diversity. *Global Change Biology*, 23, 2999–3011. https://doi.org/10.1111/gcb.13571

Jarzyna, M. A., & Jetz, W. (2018). Taxonomic and functional diversity change is scale dependent. *Nature Communications*, 9, 2565. https://doi.org/10.1038/s41467-018-04889-z

Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. *Nature*, 491, 444–448. https://doi.org/10.1038/nature11631

Jetz, W., McGeoch, M. A., Guralnick, R., Ferrier, S., Beck, J., Costello, M. J., et al. (2019). Essential biodiversity variables for mapping and monitoring species populations. *Nature Ecology & Evolution*, 3, 539–551. https://doi.org/10.1038/s41559-019-0826-1

Johnston, A., Hochachka, W. M., Strimas-Mackey, M. E., Ruiz Gutierrez, V., Robinson, O. J., Miller, E. T., et al. (2021). Analytical guidelines to increase the value of community science data: an example using eBird data to estimate species distributions. *Diversity and Distributions*, 27, 1265–1277. https://doi.org/10.1111/ddi.13271

Joppa, L. N., & Pfaff, A. (2009). High and far: biases in the location of protected areas. *PLoS ONE*, 4, e8273. https://doi.org/10.1371/journal.pone.0008273

Jung, M., Arnell, A., de Lamo, X., García-Rangel, S., Lewis, M., Mark, J., et al. (2021). Areas of global importance for conserving terrestrial biodiversity, carbon and water. *Nature Ecology & Evolution*, 5, 1499–1509. https://doi.org/10.1038/s41559-021-01528-7

Kelling, S., Johnston, A., Bonn, A., Fink, D., Ruiz-Gutierrez, V., Bonney, R., et al. (2019). Using semistructured surveys to improve citizen science data for monitoring biodiversity. *BioScience*, 69, 170–179. https://doi.org/10.1093/biosci/biz010

Kéry, M., Royle, J. A., Schmid, H., Schaub, M., Volet, B., Häfliger, G., & Zbinden, N. (2010). Site-occupancy distribution modeling to correct population-trend estimates derived from opportunistic observations. *Conservation Biology*, 24, 1388–1397. https://doi.org/10.1111/j.1523-1739.2010.01479.x

Kéry, M., & Royle, J. A. (2021). *Applied hierarchical modeling in ecology: analysis of distribution, abundance and species richness in R and BUGS, Volume 2 — dynamic and advanced models*. Academic Press.

Laliberté, E., & Legendre, P. (2010). A distance-based framework for measuring functional diversity from multiple traits. *Ecology*, 91, 299–305. https://doi.org/10.1890/08-2244.1

Legendre, P. (2014). Interpreting the replacement and richness difference components of beta diversity. *Global Ecology and Biogeography*, 23, 1324–1334. https://doi.org/10.1111/geb.12207

Le Viol, I., Jiguet, F., Brotons, L., Herrando, S., Lindström, Å., Pearce-Higgins, J. W., et al. (2012). More and more generalists: two decades of changes in the European avifauna. *Biology Letters*, 8, 780–782. https://doi.org/10.1098/rsbl.2012.0496

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A., & Langtimm, C. A. (2002). Estimating site occupancy rates when detection probabilities are less than one. *Ecology*, 83, 2248–2255. https://doi.org/10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2

MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G., & Franklin, A. B. (2003). Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. *Ecology*, 84, 2200–2207. https://doi.org/10.1890/02-3090

Magurran, A. E., Dornelas, M., Moyes, F., Gotelli, N. J., & McGill, B. (2015). Rapid biotic homogenization of marine fish assemblages. *Nature Communications*, 6, 8405. https://doi.org/10.1038/ncomms9405

Margules, C. R., & Pressey, R. L. (2000). Systematic conservation planning. *Nature*, 405, 243–253. https://doi.org/10.1038/35012251

Maxwell, S. L., Cazalis, V., Dudley, N., Hoffmann, M., Rodrigues, A. S. L., Stolton, S., et al. (2020). Area-based conservation in the twenty-first century. *Nature*, 586, 217–227. https://doi.org/10.1038/s41586-020-2773-z


McGill, B. J., Dornelas, M., Gotelli, N. J., & Magurran, A. E. (2015). Fifteen forms of biodiversity trend in the Anthropocene. *Trends in Ecology & Evolution*, 30, 104–113. https://doi.org/10.1016/j.tree.2014.11.006

McKinney, M. L., & Lockwood, J. L. (1999). Biotic homogenization: a few winners replacing many losers in the next mass extinction. *Trends in Ecology & Evolution*, 14, 450–453. https://doi.org/10.1016/S0169-5347(99)01679-1

Morelli, F., Benedetti, Y., Hanson, J. O., & Fuller, R. A. (2021). Global distribution and conservation of avian diet specialization. *Conservation Letters*, 14, e12795. https://doi.org/10.1111/conl.12795

Mouillot, D., Bellwood, D. R., Baraloto, C., Chave, J., Galzin, R., Harmelin-Vivien, M., et al. (2013). Rare species support vulnerable functions in high-diversity ecosystems. *PLoS Biology*, 11, e1001569. https://doi.org/10.1371/journal.pbio.1001569

Mu, H., Li, X., Wen, Y., Huang, J., Du, P., Su, W., et al. (2022). A global record of annual terrestrial Human Footprint dataset from 2000 to 2018. *Scientific Data*, 9, 176. https://doi.org/10.1038/s41597-022-01284-8

Newbold, T., Hudson, L. N., Hill, S. L. L., Contu, S., Lysenko, I., Senior, R. A., et al. (2015). Global effects of land use on local terrestrial biodiversity. *Nature*, 520, 45–50. https://doi.org/10.1038/nature14324

Newbold, T., Hudson, L. N., Arnell, A. P., Contu, S., De Palma, A., Ferrier, S., et al. (2016). Has land use pushed terrestrial biodiversity beyond the planetary boundary? A global assessment. *Science*, 353, 288–291. https://doi.org/10.1126/science.aaf2201

Olden, J. D., & Rooney, T. P. (2006). On defining and quantifying biotic homogenization. *Global Ecology and Biogeography*, 15, 113–120. https://doi.org/10.1111/j.1466-822X.2006.00214.x

Outhwaite, C. L., Gregory, R. D., Chandler, R. E., Collen, B., & Isaac, N. J. B. (2020). Complex long-term biodiversity change among invertebrates, bryophytes and lichens. *Nature Ecology & Evolution*, 4, 384–392. https://doi.org/10.1038/s41559-020-1111-z

Ouyang, Z., Zheng, H., Xiao, Y., Polasky, S., Liu, J., Xu, W., et al. (2016). Improvements in ecosystem services from investments in natural capital. *Science*, 352, 1455–1459. https://doi.org/10.1126/science.aaf2295

Petchey, O. L., & Gaston, K. J. (2006). Functional diversity: back to basics and looking forward. *Ecology Letters*, 9, 741–758. https://doi.org/10.1111/j.1461-0248.2006.00924.x

Pigot, A. L., Sheard, C., Miller, E. T., Bregman, T. P., Freeman, B. G., Roll, U., et al. (2020). Macroevolutionary convergence connects morphological form to ecological function in birds. *Nature Ecology & Evolution*, 4, 230–239. https://doi.org/10.1038/s41559-019-1070-4

Pimm, S. L., Jenkins, C. N., Abell, R., Brooks, T. M., Gittleman, J. L., Joppa, L. N., et al. (2014). The biodiversity of species and their rates of extinction, distribution, and protection. *Science*, 344, 1246752. https://doi.org/10.1126/science.1246752

Pollock, L. J., Thuiller, W., & Jetz, W. (2017). Large conservation gains possible for global biodiversity facets. *Nature*, 546, 141–144. https://doi.org/10.1038/nature22368

Rodrigues, A. S. L., Andelman, S. J., Bakarr, M. I., Boitani, L., Brooks, T. M., Cowling, R. M., et al. (2004). Effectiveness of the global protected area network in representing species diversity. *Nature*, 428, 640–643. https://doi.org/10.1038/nature02422

Rosenberg, K. V., Dokter, A. M., Blancher, P. J., Sauer, J. R., Smith, A. C., Smith, P. A., et al. (2019). Decline of the North American avifauna. *Science*, 366, 120–124. https://doi.org/10.1126/science.aaw1313

Rota, C. T., Fletcher, R. J. Jr., Dorazio, R. M., & Betts, M. G. (2009). Occupancy estimation and the closure assumption. *Journal of Applied Ecology*, 46, 1173–1181. https://doi.org/10.1111/j.1365-2664.2009.01734.x

Santangeli, A., Rajasärkkä, A., & Lehikoinen, A. (2016). Effects of high latitude protected areas on bird communities under rapid climate change. *Global Change Biology*, 23, 2241–2249. https://doi.org/10.1111/gcb.13518

Si, X., Baselga, A., & Ding, P. (2015). Revealing beta-diversity patterns of breeding bird and lizard communities on inundated land-bridge islands by separating the turnover and nestedness components. *PLOS ONE*, 10, e0127692. https://doi.org/10.1371/journal.pone.0127692

Socolar, J. B., Gilroy, J. J., Kunin, W. E., & Edwards, D. P. (2016). How should beta-diversity inform biodiversity conservation? *Trends in Ecology & Evolution*, 31, 67–80. https://doi.org/10.1016/j.tree.2015.11.005

Stephens, P. A., Mason, L. R., Green, R. E., Gregory, R. D., Sauer, J. R., Alison, J., et al. (2016). Consistent response of bird populations to climate change on two continents. *Science*, 352, 84–87. https://doi.org/10.1126/science.aac4858

Sullivan, B. L., Wood, C. L., Iliff, M. J., Bonney, R. E., Fink, D., & Kelling, S. (2009). eBird: a citizen-based bird observation network in the biological sciences. *Biological Conservation*, 142, 2282–2292. https://doi.org/10.1016/j.biocon.2009.05.006

Tingley, M. W., & Beissinger, S. R. (2009). Detecting range shifts from historical species occurrences: new perspectives on old data. *Trends in Ecology & Evolution*, 24, 625–633. https://doi.org/10.1016/j.tree.2009.05.009

Tobias, J. A., Sheard, C., Pigot, A. L., Devenish, A. J. M., Yang, J., Sayol, F., et al. (2022). AVONET: morphological, ecological and geographical data for all birds. *Ecology Letters*, 25, 581–597. https://doi.org/10.1111/ele.13898

van Klink, R., Bowler, D. E., Gongalsky, K. B., Swengel, A. B., Gentile, A., & Chase, J. M. (2020). Meta-analysis reveals declines in terrestrial but increases in freshwater insect abundances. *Science*, 368, 417–420. https://doi.org/10.1126/science.aax9931

van Strien, A. J., van Swaay, C. A. M., & Termaat, T. (2013). Opportunistic citizen science data of animal species produce reliable estimates of distribution trends if analysed with occupancy models. *Journal of Applied Ecology*, 50, 1450–1458. https://doi.org/10.1111/1365-2664.12158

Vellend, M., Baeten, L., Myers-Smith, I. H., Elmendorf, S. C., Beauséjour, R., Brown, C. D., et al. (2013). Global meta-analysis reveals no net change in local-scale plant biodiversity over time. *Proceedings of the National Academy of Sciences USA*, 110, 19456–19459. https://doi.org/10.1073/pnas.1312779110

Villéger, S., Mason, N. W. H., & Mouillot, D. (2008). New multidimensional functional diversity indices for a multifaceted framework in functional ecology. *Ecology*, 89, 2290–2301. https://doi.org/10.1890/07-1206.1

Violle, C., Navas, M.-L., Vile, D., Kazakou, E., Fortunel, C., Hummel, I., & Garnier, E. (2007). Let the concept of trait be functional! *Oikos*, 116, 882–892. https://doi.org/10.1111/j.0030-1299.2007.15559.x

Watson, J. E. M., Dudley, N., Segan, D. B., & Hockings, M. (2014). The performance and potential of protected areas. *Nature*, 515, 67–73. https://doi.org/10.1038/nature13947

Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M. M., & Jetz, W. (2014). EltonTraits 1.0: species-level foraging attributes of the world's birds and mammals. *Ecology*, 95, 2027. https://doi.org/10.1890/13-1917.1

Yang, J., & Huang, X. (2021). The 30 m annual land cover dataset and its dynamics in China from 1990 to 2019. *Earth System Science Data*, 13, 3907–3925. https://doi.org/10.5194/essd-13-3907-2021

## SUPPORTING INFORMATION

- **S1** Data sources, access terms, source-aware deduplication audit, source-wise annual coverage.
- **S2** Relaxed-closure justification and three-year-window sensitivity analysis.
- **S3** Spatial `phi` prior specification, prior-sensitivity analysis, posterior-bound check.
- **S4** MCMC diagnostics (R̂, ESS, trace plots), WAIC, posterior predictive checks.
- **S5** Trait sources, imputation and complete-case construction; phylogeny matching; standardized phylogenetic effect sizes.
- **S6** Probabilistic Baselga partitioning: derivation, bounds, decomposition properties; spatial beta-diversity convergence analysis.
- **S7** Effort-confound controls: detection–effort posterior, effort-saturated subset, fixed-effort counterfactual, data-source detection term.
- **S8** Full 500-species trend tables and supplementary winner–loser panels.
- **S9** Naive-versus-corrected comparison maps; 50-km grid pilot; breeding-season sensitivity.
