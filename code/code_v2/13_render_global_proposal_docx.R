#!/usr/bin/env Rscript
## 13_render_global_proposal_docx.R
##
## 仿照 ~/Desktop/Proposal for Newton Fellowship-0225.docx 的结构，
## 写一份"全球尺度多物种动态占域 × 鸟类群落重组"的研究 proposal。
## 输出：results_v2/proposal_global_bird_dynamic_occupancy.docx

suppressPackageStartupMessages({ library(officer) })

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
P <- ensure_v2_dirs()

H1 <- function(d, t) body_add_par(d, t, style = "heading 1")
H2 <- function(d, t) body_add_par(d, t, style = "heading 2")
PP <- function(d, t) body_add_par(d, t, style = "Normal")
BR <- function(d) body_add_par(d, "", style = "Normal")

doc <- read_docx()

doc <- doc |>
  body_add_par("Detection-corrected global reorganisation of bird communities under accelerating environmental change",
               style = "heading 1") |>
  body_add_par("A multi-source, multi-species dynamic occupancy framework for citizen-science avian data (eBird, GBIF, iNaturalist)",
               style = "heading 2") |>
  BR()

## --- 1. Context ----------------------------------------------------------
doc <- H1(doc, "1. Context")
doc <- PP(doc, paste0(
  "Citizen-science platforms have transformed our ability to document avian biodiversity at planetary scale: ",
  "eBird now hosts >1.6 billion observations, GBIF aggregates >2 billion vertebrate records, and iNaturalist ",
  "contributes a rapidly expanding stream of research-grade observations spanning previously under-sampled ",
  "regions (Sullivan et al. 2014; Heberling et al. 2021; Callaghan et al. 2024). These data are central to ",
  "the Global Biodiversity Framework's monitoring of indicator species and to recent flagship analyses of ",
  "biodiversity change (Sun et al. 2022; He et al. 2025; Lu et al. 2026). However, they share three pervasive ",
  "biases that critically distort headline trends: imperfect and heterogeneous detection probability, ",
  "explosive growth in observer effort over the past decade, and strong spatial inequality in coverage ",
  "(Callaghan et al. 2024; Johnston et al. 2023). Naive richness, occurrence-frequency, or even ",
  "abundance-based summaries therefore conflate true ecological change with observational expansion - a ",
  "problem that has been highlighted but is rarely solved at global scale (Boyd et al. 2023; He et al. 2025)."
))
doc <- PP(doc, paste0(
  "Hierarchical multi-species dynamic occupancy models (Doser et al. 2022; MacKenzie et al. 2003) offer a ",
  "principled solution by jointly estimating the latent occupancy probability and the conditional detection ",
  "probability across repeat visits, partial-pooling species-level effects within a Bayesian community. ",
  "Recent extensions integrate AR1 temporal structure, joint posterior propagation to community indices, and ",
  "Bayesian regression diagnostics, enabling rigorous trend inference where naive analyses fail (Doser et ",
  "al. 2024; Johnston et al. 2023). Yet despite the abundance of data, three gaps remain. First, no analysis ",
  "to date has fitted a coherent multi-species multi-season occupancy model to the merged eBird + GBIF + ",
  "iNaturalist corpus at a global, realm-stratified resolution. Second, existing global avian trend studies ",
  "have not propagated detection uncertainty to phylogenetic and functional diversity, leaving community ",
  "reassembly invisible. Third, biotic homogenisation - documented at national scale (Newbold et al. 2019; ",
  "He et al. 2025) - has not been quantified globally with detection correction or causally linked to ",
  "anthropogenic pressure gradients. I will close all three gaps."
))
doc <- PP(doc, paste0(
  "Building on a fully reproducible pilot for China (1,308 grids, 200 species, 4 chains, R-hat <= 1.09; ",
  "biotic homogenisation declined ~20% across five 5-year periods after detection correction; methods, code, ",
  "and figures available in this project), the global extension will provide the first detection-corrected, ",
  "multi-dimensional global assessment of how avian communities have reorganised across the citizen-science ",
  "era (2000-2024) and which environmental drivers most strongly shape their dynamics."
))
doc <- BR(doc)

## --- 2. Objectives -------------------------------------------------------
doc <- H1(doc, "2. Objectives")
doc <- PP(doc, paste0(
  "To quantify how 25 years of accelerating land-use change, urbanisation, and climate warming have ",
  "reshaped avian biodiversity globally, after explicitly controlling for detection bias and observer-effort ",
  "expansion in citizen-science data, and moving beyond species richness to examine taxonomic, ",
  "phylogenetic, and functional community reassembly across realms and biomes."
))
doc <- PP(doc, paste0(
  "I will deliver: (1) a detection-corrected global atlas of bird occupancy and richness on a 25 km equal-",
  "area grid (zoogeographic-realm stratified), with full posterior credible intervals; (2) global maps and ",
  "decompositions of temporal beta diversity (Bray-Curtis, Sorensen with Baselga turnover and nestedness) and ",
  "biotic homogenisation; (3) trait- and functional-group rules identifying global \"winner\" and \"loser\" ",
  "species under combined pressures; and (4) a fully reproducible open-source pipeline (R + Stan/cmdstanr) ",
  "extending and integrating the spOccupancy and brms ecosystems for citizen-science-scale inference."
))
doc <- BR(doc)

## --- 3. Research questions, hypotheses, predictions ----------------------
doc <- H1(doc, "3. Research questions, hypotheses and predictions")
doc <- H2(doc, "Research questions (Q)")
doc <- PP(doc, paste0(
  "Q1 Headline trends: After detection correction, how have global bird occupancy, taxonomic richness, ",
  "phylogenetic diversity (PD) and functional trait volume changed across realms and biomes between ",
  "2000-2004 and 2020-2024?"
))
doc <- PP(doc, paste0(
  "Q2 Community reorganisation: Are bird assemblages becoming more spatially homogeneous globally, and is ",
  "the dynamic dominated by species turnover (replacement) or nestedness (richness-difference) processes?"
))
doc <- PP(doc, paste0(
  "Q3 Driver attribution: How do climate warming, land-use intensity, urbanisation pressure, and ",
  "habitat heterogeneity independently and jointly explain detection-corrected richness and beta-diversity ",
  "trends, after partialling out spatial structure?"
))
doc <- H2(doc, "Hypotheses (H)")
doc <- PP(doc, paste0(
  "H1 Decoupled dimensions: Apparent gains in citizen-science-driven richness mask phylogenetic and ",
  "functional contraction in regions of high human pressure - taxonomic richness can rise while PD/MPD/",
  "trait-volume decline."
))
doc <- PP(doc, paste0(
  "H2 Global homogenisation, turnover-driven: Pairwise community dissimilarity declines across most ",
  "biogeographic realms; this is dominated by Baselga turnover (specialist replacement by generalists) ",
  "rather than nestedness (selective loss)."
))
doc <- PP(doc, paste0(
  "H3 Trait and pressure interaction: Slow-lived, narrow-niched, low-dispersal specialists are losers ",
  "globally; fast-lived, broad-niched generalists are winners; the asymmetry is amplified in tropical and ",
  "high-urbanisation grids where thermal and habitat constraints compound."
))
doc <- H2(doc, "Predictions (P)")
doc <- PP(doc, paste0(
  "P1.1 Detection-corrected richness trends will be flat or modestly positive globally, while PD and ",
  "trait-volume trends will exhibit stronger spatial heterogeneity, with negative slopes concentrated in ",
  "tropical primary forests and human-dominated landscapes."
))
doc <- PP(doc, paste0(
  "P1.2 Realm contrasts: Indomalayan and Neotropical realms show the strongest decoupling between ",
  "taxonomic and functional/phylogenetic dimensions; Palearctic temperate grids show effort-driven ",
  "richness inflation in naive but not in detection-corrected analyses."
))
doc <- PP(doc, paste0(
  "P2.1 Median pairwise Sorensen between grids declines monotonically across the five 5-year periods in ",
  ">=70% of biogeographic realms; the decline is concentrated within urbanisation hotspots and within ",
  "intensified-agriculture grids."
))
doc <- PP(doc, paste0(
  "P2.2 Baselga decomposition: turnover dominates nestedness in >=80% of grids; the residual nestedness ",
  "signal concentrates in islands and high-elevation isolates where selective loss is expected."
))
doc <- PP(doc, paste0(
  "P3.1 Across all realms, occupancy decline correlates with smaller geographic range, narrower thermal ",
  "tolerance, and lower hand-wing index (low dispersal); generalist trait suites positively respond to ",
  "human modification at moderate intensities, with diminishing returns beyond 50% built or cropland cover."
))
doc <- PP(doc, paste0(
  "P3.2 Trait-mediated species responses scale up to community signatures consistent with H2: regions where ",
  "specialists decline most strongly are those exhibiting the steepest local homogenisation."
))
doc <- BR(doc)

## --- 4. Methodology ------------------------------------------------------
doc <- H1(doc, "4. Methodology")
doc <- H2(doc, "4.1 Data and key variables")
doc <- PP(doc, paste0(
  "Biodiversity (occurrence): The full eBird Basic Dataset (~1.6B records); GBIF Aves snapshot (~2B ",
  "records, focusing on research-grade or human-observation records with valid coordinates and dates); ",
  "iNaturalist research-grade Aves observations. Records will be deduplicated against an event key (species ",
  "x date x rounded coordinates x observer) and checked against an authoritative taxonomy backbone ",
  "(IOC + Clements + GBIF Backbone, harmonised; Birdlife/HBW alias table). Pilot work in this project on a ",
  "China subset retained 7.49M unique events from 14.56M raw records (~48.6% deduplication)."
))
doc <- PP(doc, paste0(
  "Spatial-temporal design: Global 25 km equal-area grid (~76,000 land cells; Behrmann or ISEA discrete ",
  "global grid), stratified by Wallace zoogeographic realm and WWF biome. Five 5-year primary periods ",
  "(2000-2004, 2005-2009, 2010-2014, 2015-2019, 2020-2024); each calendar year within a primary period as a ",
  "secondary occasion (repeated visit) to enable detection-occupancy separation."
))
doc <- PP(doc, paste0(
  "Effort and detection covariates: visit count (unique observer-day pairs per grid-year), unique observer ",
  "count, total checklist duration (where available), checklist completeness flag (eBird), and a composite ",
  "PCA-derived effort index (PC1 of the seven effort metrics; >85% variance explained in pilot). Effort ",
  "covariates enter the detection sub-model only, not occupancy."
))
doc <- PP(doc, paste0(
  "Climate: WorldClim v2.1 BIO1-19 baselines (1970-2000) plus CRU TS V4 monthly anomalies aggregated to ",
  "25 km grids, decomposed into mean warming and heat-extreme metrics following Outhwaite et al. (2022). ",
  "Land use and pressure: ESA WorldCover annual fractions, Human Modification Index (HMI; Theobald et al. ",
  "2020), Human Footprint Index (HFI; Williams et al. 2020), urbanisation gradient from GHSL settlement ",
  "data. Topography and habitat heterogeneity: SRTM elevation mean and SD, EarthEnv land-surface texture, ",
  "MODIS NPP and NDVI."
))
doc <- PP(doc, paste0(
  "Traits and phylogeny: AVONET (Tobias et al. 2022) for body mass, hand-wing index, beak morphology, ",
  "habitat, and range size; missForest imputation for residual gaps. Phylogeny: McTavish-Miller / Cornell ",
  "Lab open avian tree of life (clootl, version 1.6, taxonomy 2025) extracted via probability-weighted Faith's ",
  "PD and abundance-weighted MPD."
))
doc <- H2(doc, "4.2 Statistical analysis")
doc <- PP(doc, paste0(
  "I will use detection-corrected, posterior-mean occupancy under the lowest-pressure environmental baseline ",
  "(primary-forest cells in pre-2010 periods, low HFI/HMI quartile) as the reference for visualising ",
  "relative change. All models are Bayesian hierarchical, fitted in parallel for each of the eight ",
  "Wallace realms to retain ecologically meaningful contrasts."
))
doc <- PP(doc, paste0(
  "Q1 (H1) - Multi-species dynamic occupancy: For each realm I will fit ",
  "spOccupancy::tMsPGOcc(occ.formula = ~ climate + land-use + topo + space, ",
  "det.formula = ~ log_events + log_observers + duration_min, ar1 = TRUE), ",
  "with 4 chains (>=4,000 sweeps each, n.thin=10, n.burn=2,000), pooling species-level coefficients ",
  "around community means. Posterior occupancy ψ is propagated to community indices (richness, Shannon, ",
  "probability-weighted Faith's PD, abundance-weighted MPD, CWM/CWSD trait moments, trait-space volume, ",
  "Rao's Q) for every retained MCMC draw, yielding mean + 95% credible intervals at the grid x period ",
  "level. Convergence is verified with Gelman-Rubin R-hat (target <= 1.05), effective sample size, and ",
  "spOccupancy::ppcOcc Freeman-Tukey Bayesian p-values. The pilot (200 species, China) achieved R-hat ",
  "median 1.007 and posterior predictive consistency, demonstrating feasibility of scaling."
))
doc <- PP(doc, paste0(
  "Q2 (H2) - Community reassembly: Per-grid temporal beta between consecutive periods is decomposed into ",
  "Bray-Curtis dissimilarity and Baselga turnover/nestedness components on probability-Sorensen rather than ",
  "binary occurrence, propagating posterior uncertainty. Spatial biotic homogenisation is quantified per ",
  "period as median pairwise Sorensen dissimilarity between each grid and its K=20 nearest geographic ",
  "neighbours; per-grid temporal slopes provide a global map of where homogenisation is strongest."
))
doc <- PP(doc, paste0(
  "Q3 (H3) - Driver and trait attribution: Per-grid 5-period trend slopes of richness, PD, and trait-",
  "volume are modelled with brms (Stan via cmdstanr backend, 4 chains, 4,000 iterations) on standardised ",
  "covariates from four driver groups (climate, topo+habitat, human pressure, pure space). Independent ",
  "variance fractions are estimated with vegan::varpart on the same matrix; DHARMa simulated residuals ",
  "diagnose dispersion and spatial autocorrelation. At the species level I will fit phylogenetically ",
  "informed brms models of occupancy trend and range shift on traits (body mass, hand-wing index, ",
  "longevity, range size, thermal-niche breadth), with phylogenetic correlation structure derived from ",
  "the clootl tree, identifying global \"winner\" and \"loser\" trait syndromes."
))
doc <- BR(doc)

## --- 5. Milestones -------------------------------------------------------
doc <- H1(doc, "5. Milestones and timescales (24 months)")

mile_tbl <- data.frame(
  Period = c("Months 1-4", "Months 5-10", "Months 11-16",
              "Months 17-20", "Months 21-24"),
  `Key activities` = c(
    "Data harmonisation: eBird/GBIF/iNaturalist taxonomic and event-level deduplication; global 25 km grid construction; effort-covariate compilation; pilot fits of tMsPGOcc on Palearctic and Nearctic realms",
    "Realm-by-realm production fits (Indomalayan, Afrotropical, Neotropical, Australasian, Oceanian, Antarctic); convergence diagnostics and PPC; Paper 1 drafting (global detection-corrected occupancy and richness atlas)",
    "Community reassembly module: posterior propagation to PD/trait diversity; Baselga turnover/nestedness decomposition; spatial homogenisation maps; pre-registered analysis plan on OSF",
    "Trait/functional-group attribution: brms phylogenetic regressions; varpart decomposition of driver groups; cross-realm comparison; Paper 2 drafting (community reassembly + trait-mediated vulnerability)",
    "Final synthesis and dissemination: open-source pipeline release; policy brief (CBD GBF / IPBES alignment); stakeholder workshops; post-fellowship roadmap"
  ),
  Deliverables = c(
    "Reproducible pipeline (GitHub); OSF pre-registration; preliminary occupancy figures",
    "Paper 1 submitted (global detection-corrected avian occupancy atlas); conference abstract (BES/ICCB)",
    "Reproducible beta-diversity workflow; key figures for Paper 2",
    "Paper 2 submitted (community reassembly and trait-mediated vulnerability)",
    "Open-science package (Zenodo + GitHub); policy brief; post-fellowship roadmap"
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
ft <- flextable::flextable(mile_tbl)
ft <- flextable::theme_vanilla(ft)
ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
doc <- flextable::body_add_flextable(doc, ft)
doc <- BR(doc)

## --- 6. Challenges, feasibility, risks -----------------------------------
doc <- H1(doc, "6. Challenges, feasibility and risks")
doc <- PP(doc, paste0(
  "Spatial coverage inequality: Citizen-science records are skewed toward Europe, North America, and East ",
  "Asia. I will report all results as realm-stratified, conduct sensitivity analyses with effort-stratified ",
  "subsamples, and treat under-sampled grids with explicit posterior caveats. Where realm-level coverage is ",
  "insufficient I will fall back to coarser biome aggregation."
))
doc <- PP(doc, paste0(
  "Computational scale: Multi-species dynamic occupancy on global 25 km grids with ~500-1,000 candidate ",
  "species per realm exceeds memory budgets of standard workstations. I will (a) fit each realm separately, ",
  "(b) thin posterior draws to <=200 retained samples for downstream propagation, (c) use cmdstanr ",
  "multi-threading where applicable, and (d) leverage the pilot infrastructure (already verified on the ",
  "China 200-species, 1,308-grid case) for predictable scaling."
))
doc <- PP(doc, paste0(
  "Detection-occupancy identifiability: tMsPGOcc requires sufficient repeat visits per grid-year. I will ",
  "exclude grid-years with <2 visits and apply weakly informative priors on detection covariates to ",
  "stabilise low-effort estimates. Posterior predictive checks (Freeman-Tukey discrepancy, Bayesian ",
  "p-value) will be reported per realm."
))
doc <- PP(doc, paste0(
  "Observational vs causal inference: Citizen-science records are observational; I will explicitly frame ",
  "results as detection-corrected associations, complement them with variance partitioning to attribute ",
  "independent driver fractions, and discuss the limits to causal claim. Spatial residuals will be ",
  "diagnosed via DHARMa and Moran's I; if residual autocorrelation is non-trivial, I will add ",
  "Gaussian-process or CAR/SAR random fields in brms as a robustness check."
))
doc <- BR(doc)

## --- 7. Expected outputs --------------------------------------------------
doc <- H1(doc, "7. Expected outputs")
doc <- PP(doc, paste0(
  "Two first-authored papers targeting Nature-family / Global Change Biology / Ecology Letters: (i) the ",
  "first global detection-corrected dynamic-occupancy atlas of bird community structure, and (ii) global ",
  "biotic homogenisation and trait-mediated reassembly across the citizen-science era. Additional outputs: ",
  "a pre-registered analysis plan (OSF); fully open code and derived datasets (Zenodo + GitHub); a ",
  "policy-facing brief aligned with the Global Biodiversity Framework and IPBES indicators; and ",
  "presentations at BES, ICCB, and the IBS World Biogeography Symposium."
))
doc <- BR(doc)

## --- 8. References --------------------------------------------------------
doc <- H1(doc, "8. References (selected)")
refs <- c(
  "Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. Glob. Ecol. Biogeogr. 19, 134-143.",
  "Boyd, R. J. et al. (2023). ROBITT: a tool for assessing the risk of bias in studies using uncoordinated biological records. Methods Ecol. Evol.",
  "Callaghan, C. T. et al. (2024). The role of citizen science in biodiversity science. Trends Ecol. Evol.",
  "Doser, J. W., Finley, A. O. & Banerjee, S. (2022). spOccupancy: An R package for single-species, multi-species, and integrated occupancy models. Methods Ecol. Evol. 13, 1670-1678.",
  "Doser, J. W. et al. (2024). spOccupancy v2: hierarchical Bayesian dynamic occupancy with AR1 temporal structure.",
  "He, J. et al. (2025). Decade-long bird trends in China: stable species richness but increasing biotic homogenisation. Glob. Change Biol.",
  "Heberling, J. M. et al. (2021). Data integration enables global biodiversity synthesis. PNAS 118, e2018093118.",
  "Johnston, A. et al. (2023). Best practices for making reliable inferences from citizen science data. Methods Ecol. Evol.",
  "Lu, X. et al. (2026). Multidecadal legacy of uneven urbanisation on divergent prospects for bird biodiversity. Nat. Cities 3, 176-188.",
  "MacKenzie, D. I. et al. (2003). Estimating site occupancy, colonisation, and local extinction when a species is detected imperfectly. Ecology 84, 2200-2207.",
  "Newbold, T. et al. (2019). Climate and land-use change homogenise terrestrial biodiversity. Emerg. Top. Life Sci. 3, 207-219.",
  "Outhwaite, C. L., McCann, P. & Newbold, T. (2022). Agriculture and climate change are reshaping insect biodiversity worldwide. Nature 605, 97-102.",
  "Sullivan, B. L. et al. (2014). The eBird enterprise: an integrated approach to development and application of citizen science. Biol. Conserv. 169, 31-40.",
  "Sun, B. et al. (2022). Urbanisation affects spatial variation and species similarity of bird diversity distribution. Sci. Adv. 8, eabq9212.",
  "Theobald, D. M. et al. (2020). Earth transformed: detailed mapping of global human modification. Earth Syst. Sci. Data 12, 1953-1972.",
  "Tobias, J. A. et al. (2022). AVONET: morphological, ecological and geographical data for all birds. Ecol. Lett. 25, 581-597.",
  "Williams, J. J. & Newbold, T. (2020). Local climatic changes affect biodiversity responses to land use. Divers. Distrib. 26, 76-92."
)
for (r in refs) doc <- PP(doc, r)

out_path <- v2_file("results", "proposal_global_bird_dynamic_occupancy", "docx")
print(doc, target = out_path)
message(sprintf("[stage-13] global proposal: %s", out_path))
