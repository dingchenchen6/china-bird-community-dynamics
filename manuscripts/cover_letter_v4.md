# Cover letter

**Subject:** Submission of "Detection-corrected community dynamics of Chinese birds (2000–2024)"

Dear Editor,

We are pleased to submit our manuscript "Detection-corrected community dynamics of Chinese birds, 2000–2024" for consideration at **Global Change Biology** / **Nature Communications** / **Methods in Ecology and Evolution** (please adjust based on target journal).

**Why this paper, why now.** Citizen-science platforms now provide multi-decadal, country-scale biodiversity data unattainable from any single research programme, but the community is increasingly aware that uncorrected occurrence frequencies can manufacture spurious trends purely from sampling expansion. Our study is the first to integrate the China Birdwatching Records Center (the largest national platform) with eBird/GBIF China under a transparent source-aware deduplication, and to fit a **spatial multi-species dynamic occupancy model** at country scale across 25 years. We propagate posterior occupancy uncertainty consistently into 10 community metrics (taxonomic, functional, **probability-weighted Faith's PD**), into Baselga's turnover-vs-nestedness decomposition, and into per-species trends and trait-mediated predictors — yielding the first fully Bayesian, uncertainty-propagated portrait of how Chinese avian communities have reshaped under intertwined climate, land-use and human-pressure change.

**Three novel contributions.**

1. **Methodological**: We demonstrate that the conventional naive-richness approach flips the trend direction in {nvc_flip_rate}% of grids relative to detection-corrected richness, underscoring how citizen-science conclusions can be systematically biased without occupancy correction.
2. **Ecological**: We quantify, at country scale, that Chinese bird community β-diversity is **turnover-dominated** rather than nestedness-driven, with the dominance varying spatially across climatic gradients and human-pressure zones.
3. **Trait-based prediction**: By integrating AVONET, EltonTraits and IUCN habitat breadth, we identify a coherent expander phenotype (high HWI, broad habitat breadth, low diet specialization), reinforced by phylogenetic constraints (LOO compares M0 vs M1).

**Robustness.** We provide four sensitivity analyses (3-yr vs 5-yr window; breeding-season vs year-round; 100-km vs 50-km grid; eps threshold from 1e-12 to 0.10) demonstrating that conclusions are robust to design choices. All code and derived data are deposited on Zenodo with DOI; raw occurrence records are shared under provider agreements.

**Open science.** Full reproducibility is supported: `code_v4/` is modular (data merge → survey history → environment → stMsPGOcc → diagnostics → multi-diversity → driver regression → trait regression → publication figures), each stage seeded with `set_seeds()`, and large posteriors are stored via `qs` to minimise reviewer footprint.

We believe this work meets the readership of [target journal] because it both delivers a continent-scale empirical finding on biotic homogenization under global change, and supplies an open, well-tested methodological pipeline that other regions can readily adopt.

We confirm that this manuscript has not been published or submitted elsewhere, that all authors have approved the submission, and that we have no competing interests to declare. Five suggested reviewers are listed at the end of this letter.

Thank you for your consideration.

Sincerely,

**Chenchen Ding** (corresponding author)
School of Life Sciences / Institute of Ecology, Peking University
chenchending1992 [at] gmail.com

---

### Suggested reviewers
1. Jeffrey W. Doser (Michigan State University) — author of spOccupancy
2. Mark A. Tucker (Senckenberg) — human-pressure × wildlife movement expertise
3. Andrés Baselga (USC Spain) — β-diversity decomposition expert
4. Frédéric Jiguet (MNHN, Paris) — citizen-science bird community trends in Europe
5. Yiming Hu (Sun Yat-sen University) — Chinese bird biogeography expert
