# Hidden Homogenization in a Warming Megadiverse Nation

## A Detection-Corrected Early-Warning System for Chinese Bird Community Change

**Proposal type:** High-impact research proposal / fellowship / foundation grant / flagship collaborative project  
**Prepared from:** local China bird dynamic occupancy analysis, v3 result tables, v4 audit materials  
**Date:** 2026-06-02  
**Lead concept:** Turn a national citizen-science bird record archive into a detection-corrected biodiversity intelligence system that separates real ecological change from observation expansion.

## 1. Executive Summary

China now has enough citizen-science bird data to ask a question that was impossible a decade ago: are bird communities truly recovering, reorganizing, or merely becoming easier to observe? Raw occurrence records suggest strong increases in bird richness across many regions, but citizen-science participation, digital reporting, access, platform protocols, and observer effort have changed dramatically since 2000. Without explicit detection correction, biodiversity gain can be confused with monitoring expansion.

This project proposes a rigorous, scalable framework for national bird community assessment based on dynamic multispecies occupancy models. We will integrate China Birdwatching Records Center data, eBird/GBIF records, environmental change layers, functional traits, and phylogenetic information to quantify 25 years of bird community change across China. The core innovation is not simply mapping species trends. It is propagating detection-corrected occupancy uncertainty into taxonomic, functional, phylogenetic, and beta-diversity metrics, then identifying where apparent biodiversity gains mask functional homogenization.

Our proof-of-concept is already compelling. A 200-species spatial dynamic occupancy analysis shows mean corrected richness rising from 79.15 to 100.83 species per grid from 2000-2004 to 2020-2024. A 500-species temporal extension shows the same pattern at broader taxonomic scale, with corrected richness increasing from 155.04 to 195.55 species per grid. Yet functional trait volume and Rao's Q decline slightly, and beta diversity shifts sharply from turnover to nestedness in the latest interval. These findings reveal a powerful conservation message: China's bird communities may be gaining species while losing distinctiveness.

The proposed project will transform this proof-of-concept into a submission-ready scientific product, an open reproducible analysis system, and a national monitoring template for detection-corrected biodiversity intelligence.

## 2. The Big Idea

### From observation growth to ecological truth

Citizen science is often treated as a solution to biodiversity data scarcity. In reality, it creates a new scientific challenge: the observation process itself changes through time. More observers, better mobile apps, improved identification skills, and increasing travel can all produce apparent biodiversity gains.

The central idea of this proposal is to build an inference system that asks, for every grid cell, period, and species:

1. What was probably present?
2. What was probably detected?
3. How uncertain are those estimates?
4. What happens to community-level diversity after that uncertainty is propagated?

This turns citizen-science data from a record of where people went birding into an ecological instrument for measuring where bird communities are changing.

### The conceptual advance

Most large-scale biodiversity assessments report taxonomic richness, occupancy trends, or raw species lists. This project integrates four layers:

- **Detection correction:** dynamic multispecies occupancy models separate latent occupancy from observation effort.
- **Multidimensional diversity:** posterior occupancy probabilities are propagated into richness, Shannon diversity, functional structure, phylogenetic structure, and beta diversity.
- **Community reorganization:** probability-weighted Baselga decomposition separates turnover from nestedness over time.
- **Mechanistic interpretation:** species traits and environmental gradients identify likely winners, losers, and homogenization pathways.

The result is a new class of national biodiversity early-warning system: one that can distinguish "more species observed" from "communities becoming more similar."

## 3. Why This Matters Now

China is undergoing rapid ecological change: climate warming, wetland restoration, urban expansion, agricultural intensification, reforestation, infrastructure growth, and protected-area reform all interact across a vast environmental gradient. Birds are excellent sentinels of these changes because they are mobile, observable, ecologically diverse, and culturally visible.

At the same time, China has experienced explosive growth in citizen-science bird observation. This creates both opportunity and danger. The opportunity is unprecedented national coverage. The danger is that observation expansion can masquerade as biological recovery.

This project is timely because it addresses both sides simultaneously. It leverages the national growth of bird data while explicitly correcting for the biases introduced by that same growth.

## 4. Proof-of-Concept Evidence

### 4.1 Spatial core: 200-species dynamic occupancy analysis

The primary spatial analysis covers 200 bird species across 1,247 100-km grid cells and five 5-year periods.

Key proof-of-concept findings:

- Mean occupancy-corrected richness rises from 79.15 to 100.83 species per grid.
- Shannon diversity rises from 4.74 to 4.88.
- Inverse Simpson diversity rises from 106.37 to 126.64.
- Functional trait volume declines from 1.383 to 1.365.
- Rao's Q declines from 1.705 to 1.680.
- Species trends are expansion dominated: 136 expanding, 62 stable, 2 contracting.
- Naive occurrence trends flip direction after occupancy correction in 37 of 1,247 grids.
- Beta-diversity turnover proportion declines from 65.9% in P1-P2 to 11.3% in P4-P5.

Interpretation: corrected richness increases are real enough to survive detection correction, but the diversity added is functionally conservative and increasingly nested.

### 4.2 Breadth test: 500-species temporal extension

The 500-species extension covers the same grid and period structure using a temporal multispecies dynamic occupancy model.

Key proof-of-concept findings:

- Mean corrected richness rises from 155.04 to 195.55 species per grid.
- Shannon rises from 5.507 to 5.645.
- Inverse Simpson rises from 218.69 to 260.33.
- Functional trait volume declines from 1.581 to 1.565.
- Rao's Q declines from 1.732 to 1.711.
- Species trends: 258 expanding, 227 stable, 15 contracting.
- Naive-vs-corrected trend direction flips in 25 of 1,247 grids.
- Turnover proportion declines from 75.7% in P1-P2 to 19.5% in P4-P5.

Interpretation: the expansion-with-homogenization signal generalizes beyond the 200-species spatial subset.

### 4.3 Mechanistic clues

Complete-case trait and environmental analyses currently cover 200 species:

- Habitat breadth is positively associated with occupancy trends.
- Diet specialization is weakly negatively associated with occupancy trends.
- Cropland, elevation, temperature seasonality, wettest-month precipitation, human footprint, and built land cover are positively associated with species trend slopes.
- Random forest importance suggests that spatial baseline structure dominates, followed by land-use and climate-change variables.

Interpretation: broad-habitat and human-modified-landscape-tolerant species appear to be central drivers of national taxonomic gain.

## 5. Objectives and Hypotheses

### Objective 1. Produce a submission-ready national assessment of detection-corrected bird community dynamics

**Hypothesis 1:** Occupancy-corrected richness has increased across much of China, but the increase is smaller and spatially more structured than naive occurrence trends suggest.

**Deliverables:** verified 4-chain spatial model; convergence diagnostics; corrected-richness, Shannon, functional, phylogenetic and beta-diversity maps; manuscript-ready figures.

### Objective 2. Test whether taxonomic gains correspond to functional and phylogenetic differentiation

**Hypothesis 2:** Taxonomic diversity increases are accompanied by weak or negative changes in functional trait volume, indicating homogenization rather than multidimensional recovery.

**Deliverables:** functional diversity trajectories; McTavish/Jetz phylogenetic diversity validation; functional-phylogenetic mismatch maps.

### Objective 3. Identify species winners, losers, and trait filters

**Hypothesis 3:** Species with broader habitat breadth, lower diet specialization, greater mobility, and greater tolerance of human-modified landscapes show more positive occupancy trends.

**Deliverables:** species trend atlas; trait-trend models; winner-loser forest plots; complete-case and imputation-sensitive trait analyses.

### Objective 4. Separate observation effects from ecological effects

**Hypothesis 4:** Occupancy correction will preserve the national richness-increase direction but substantially reduce or spatially reallocate trend magnitudes, especially in regions with rapid observer growth.

**Deliverables:** naive-vs-corrected trend maps; source-bias diagnostics; platform-specific detection sensitivity.

### Objective 5. Build a reusable national biodiversity intelligence pipeline

**Hypothesis 5:** A standardized detection-corrected pipeline can convert citizen-science records into a repeatable monitoring system for annual or 5-year biodiversity reporting.

**Deliverables:** documented R pipeline; data manifest; reproducibility archive; policy-facing dashboard prototype.

## 6. Work Plan

### Work Package 1. Evidence consolidation and reproducibility hardening

**Goal:** turn the current local proof-of-concept into an auditable analysis archive.

Tasks:

- Sync or regenerate full 200sp and 500sp model fit objects.
- Store model fit checksums and code commit hashes.
- Generate `sessionInfo()` and dependency lock files.
- Create a source-to-output manifest mapping every figure/table to scripts and inputs.
- Confirm that v3 and v4 outputs are not mixed in final reporting.

Success criteria:

- Every manuscript number can be traced to a CSV, script, model object, and run label.
- Full model objects are available locally or in a secure archive.

### Work Package 2. Final spatial occupancy modelling

**Goal:** produce defensible spatial inference.

Tasks:

- Rerun or verify 200sp `stMsPGOcc` with four chains.
- Report R-hat and ESS by parameter group.
- Run posterior predictive checks.
- Test sensitivity to `phi` prior/range scale.
- Test NNGP neighbour sensitivity in pilot runs.
- Add source-specific detection structure where identifiable.

Success criteria:

- Max R-hat <= 1.05 for critical parameters or transparent explanation of exceptions.
- PPC Bayesian p-values are not systematically extreme.
- Spatial prior and NNGP settings are justified with sensitivity evidence.

### Work Package 3. 500-species breadth extension

**Goal:** elevate the 500-species result from "additional table set" to a strong generality analysis.

Tasks:

- Verify 500sp temporal model metadata and chain status.
- Generate 500sp figures for richness trajectory, species trend classes, winner-loser species, beta turnover/nestedness, and naive-vs-corrected comparison.
- Expand trait joins beyond current 200 complete cases where feasible.
- Clearly label the 500sp model as temporal/non-spatial unless a spatial version is successfully run.

Success criteria:

- 500sp results support manuscript generality claims with visual evidence, not only tables.
- Trait mechanism claims report the correct complete-case n.

### Work Package 4. Multidimensional diversity and homogenization analysis

**Goal:** make "expansion without differentiation" the central ecological contribution.

Tasks:

- Repair or replace grey PD panels using verified McTavish phylogenetic metrics.
- Produce clean maps for corrected richness, Shannon, functional trait volume, Rao's Q, FEve/FDiv, PD, and MPD.
- Generate probability-weighted Baselga maps and proportions.
- Include a mathematical supplement for probability-weighted beta partitioning.
- Avoid emphasizing synchrony until the boundary-value issue is resolved.

Success criteria:

- Main figures show taxonomic gain, functional contraction, and nestedness shift in one coherent visual narrative.

### Work Package 5. Mechanisms, traits, and conservation interpretation

**Goal:** connect statistical patterns to ecological meaning and conservation action.

Tasks:

- Fit trait-trend models with and without phylogenetic random effects.
- Compare complete-case, imputed, and expanded-trait models.
- Quantify regional priority areas: grids with corrected richness decline, high functional loss, high nestedness, or high uncertainty.
- Map observer-source bias and propose monitoring network expansion priorities.
- Produce a policy-facing brief in Chinese and English.

Success criteria:

- Manuscript claims are supported by trait/environment results without causal overreach.
- Conservation outputs identify actionable regions and monitoring gaps.

## 7. Innovation

### Statistical innovation

- Integrates dynamic multispecies occupancy modelling with posterior propagation into multidimensional diversity.
- Extends Baselga decomposition to probability-weighted occupancy inputs.
- Compares naive and detection-corrected trends at national scale.
- Combines spatial inference with broad 500-species temporal validation.

### Ecological innovation

- Reframes apparent biodiversity gain as a possible homogenization process.
- Identifies broad-habitat species as potential winners in a rapidly changing landscape.
- Links taxonomic, functional, and phylogenetic evidence in a single national assessment.

### Monitoring innovation

- Creates a reproducible template for detection-corrected biodiversity reporting.
- Converts citizen-science observation growth from a bias into a modelled component of inference.
- Provides a scalable framework for future annual or 5-year national biodiversity accounts.

## 8. Expected Outputs

### Scientific outputs

- One flagship paper targeting Nature Ecology & Evolution, Nature Communications, Global Change Biology, Ecology Letters, or Science Advances.
- One methods or data paper on probability-weighted multidimensional diversity propagation from occupancy posteriors.
- One open workflow repository with reproducible code, manifests, and figure scripts.

### Data outputs

- Detection-corrected occupancy summaries by species, grid, and period.
- Multidimensional diversity tables with credible intervals.
- Species winner-loser atlas.
- Naive-vs-corrected trend difference maps.
- Sampling and source-bias audit tables.

### Conservation outputs

- Priority maps for monitoring expansion.
- Priority maps for functional homogenization risk.
- Policy brief for national biodiversity monitoring and citizen-science data standards.

## 9. Timeline

### Months 1-2: Reproducibility and diagnostics

- Sync model objects.
- Verify v3/v4 provenance.
- Rerun 4-chain spatial model if required.
- Produce MCMC/PPC diagnostic package.

### Months 3-4: 500sp and sensitivity analyses

- Verify 500sp temporal extension.
- Generate full 500sp figure set.
- Run phi prior, NNGP neighbour, breeding-season, source detection, 3-year window, and grid-size sensitivity analyses.

### Months 5-6: Final diversity and mechanism analyses

- Repair phylogenetic metrics.
- Finalize functional/phylogenetic maps.
- Fit trait and driver models.
- Build conservation priority layers.

### Months 7-8: Writing and submission

- Finalize manuscript.
- Prepare supplementary information.
- Prepare cover letter and graphical abstract.
- Submit to target journal.

### Months 9-12: Extension and dissemination

- Build policy brief.
- Archive data/code.
- Prepare second methods/data paper.
- Present results to collaborator and conservation audiences.

## 10. Risk Management

### Risk 1. Full spatial model remains computationally expensive

Mitigation:

- Use the existing 200sp spatial model as the main spatial core.
- Use 500sp temporal model as breadth extension.
- Run spatial sensitivity on representative subsets.
- Store thinned posterior occupancy arrays separately from full model objects.

### Risk 2. Multi-chain diagnostics reveal convergence issues

Mitigation:

- Increase batch/burn-in for problematic species or parameter groups.
- Simplify AR1 or latent factor structure if not identifiable.
- Use posterior predictive checks and sensitivity analyses to identify robust conclusions.

### Risk 3. Source-specific detection effects are confounded with space

Mitigation:

- Model source where identifiable.
- Add source-proportion sensitivity in driver models.
- Map source dominance and restrict claims where inference is source-confounded.

### Risk 4. Phylogenetic matching remains incomplete

Mitigation:

- Prioritize McTavish/clootl matched tree outputs.
- Report phylogenetic metrics only for verified matched species.
- Keep functional/taxonomic results as the primary inference if phylogenetic coverage is insufficient.

### Risk 5. Trait data remain complete for only 200 species

Mitigation:

- Treat mechanism analysis as complete-case unless expanded.
- Use multiple imputation and missingness diagnostics.
- Separate trend generality from trait mechanism claims.

## 11. Why This Proposal Is Fundable

This project combines scale, urgency, methodological rigour, and policy relevance.

- **Scale:** 25 years, national extent, hundreds of species, multiple biodiversity dimensions.
- **Urgency:** China is undergoing rapid environmental change, and biodiversity monitoring must distinguish ecological change from observation growth.
- **Innovation:** detection-corrected multidimensional diversity propagation is still rare at national scale.
- **Feasibility:** a local proof-of-concept already exists with coherent 200sp and 500sp results.
- **Impact:** outputs support top-journal publication, open science infrastructure, and applied conservation monitoring.

The proposal is especially strong because it does not promise a vague future dataset. It starts from an existing working pipeline and turns it into a rigorous, auditable, publishable national assessment.

## 12. Flagship Narrative for Reviewers

> China's bird communities appear richer today than two decades ago. But is this ecological recovery, observation growth, or a subtler reorganization? By correcting imperfect detection across millions of citizen-science records and propagating uncertainty into taxonomic, functional and phylogenetic diversity, we show that richness gains can coexist with functional contraction and nested community assembly. This project will build the first detection-corrected early-warning system for Chinese bird community homogenization and create a reusable framework for biodiversity intelligence in the citizen-science era.

## 13. Target Journals and Positioning

### Primary targets

- **Nature Ecology & Evolution:** best if the final narrative emphasizes hidden homogenization, citizen-science correction, and broad ecological theory.
- **Nature Communications:** best if the full national scale, data integration, and multidimensional evidence are emphasized.
- **Global Change Biology:** best if climate/land-use change and national-scale community reorganization are emphasized.
- **Ecology Letters:** best if the conceptual novelty of expansion without differentiation and probability-weighted diversity propagation is sharpened.
- **Science Advances:** best if the manuscript is framed around national biodiversity intelligence and policy relevance.

### Positioning statement

This is not a regional bird trend paper. It is a demonstration that citizen-science biodiversity gains can conceal functional homogenization unless detectability is modelled and multidimensional community structure is propagated from posterior occupancy.

## 14. Immediate Next Steps

1. Sync or regenerate full model fit objects and logs.
2. Produce verified 4-chain convergence tables for 200sp and, if available, 500sp.
3. Repair PD figures and switch to verified McTavish outputs.
4. Create 500sp figure suite beyond RF importance.
5. Add source-specific detection audit.
6. Rewrite final manuscript figures around the "expansion without differentiation" narrative.
7. Prepare a 2-page graphical abstract and a 1-page cover letter.

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
