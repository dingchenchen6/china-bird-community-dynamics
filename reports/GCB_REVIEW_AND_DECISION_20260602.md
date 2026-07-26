# Peer-Review Simulation & Editorial Decision — *Global Change Biology*

> Manuscript reviewed: `MANUSCRIPT_NEE_v5_20260602.md` (the current best draft)
> Target journal (per author request): **Global Change Biology**
> Mode: `full` (5 reviewers + editorial synthesis). READ-ONLY — this report does not modify the manuscript; the GCB-revised manuscript is delivered separately as `MANUSCRIPT_GCB_v6_20260602.md`.
> Date: 2026-06-02

---

## Phase 0 — Field analysis & reviewer configuration

**Field analysis.** Primary discipline: global change ecology / macroecology. Secondary: Bayesian quantitative ecology (hierarchical occupancy modelling) + citizen-science data science. Paradigm: observational, model-based inference with posterior propagation. Methodology type: spatial multispecies dynamic occupancy + multidimensional biodiversity metrics + driver attribution. Target tier: top specialist (GCB, IF ~11) — a strong fit because GCB centres on biodiversity responses to climate, land-use and human-pressure change at large scales, and welcomes detailed methods and mechanism. Paper maturity: scientifically advanced, but evidence reporting (single-chain MCMC, missing CIs) and global-change framing are not yet at submission grade.

**Reviewer configuration card (GCB-calibrated). Adjustable on request.**

| # | Persona | Identity | Particular focus |
|---|---|---|---|
| EIC | Editor-in-Chief | GCB handling editor, macroecologist working on large-scale biodiversity change | Global-change relevance, conceptual advance beyond a regional case study, fit to GCB readership |
| R1 | Methodology | Bayesian hierarchical occupancy modeller (spOccupancy/JAGS/NIMBLE user) | Convergence, identifiability of detection vs occupancy, prior scaling, probabilistic beta-diversity derivation, trend statistics |
| R2 | Domain | Avian community ecologist working on biotic homogenization and functional diversity | Homogenization definition/evidence, functional-metric robustness, winners–losers, literature dialogue |
| R3 | Perspective | Biodiversity-monitoring / EBV and citizen-science science scientist | Generalizability to global monitoring, effort confound at policy scale, reproducibility/archiving |
| DA | Devil's Advocate | Skeptical quantitative ecologist | Does occupancy correction truly remove the effort artefact? Is "functional homogenization" statistically real? |

---

## Phase 1 — Independent reviews

### Reviewer EIC — Editor-in-Chief assessment

**Summary judgement.** The manuscript addresses a question of clear importance to GCB: whether continental, citizen-science-based evidence of biodiversity change reflects genuine ecological reorganization or the growth of observation itself. The "expansion without differentiation" finding is a genuine conceptual contribution and is well-matched to GCB's remit on biodiversity responses to global change. However, in its current form the paper is framed primarily as a *methods-and-region* study (occupancy correction applied to Chinese birds) rather than as a *global-change* study, which is what GCB readers expect.

**Strengths.** (1) 25-year, national-scale, detection-corrected, multidimensional analysis is ambitious and timely. (2) The taxonomic-vs-functional divergence is a real, interpretable signal. (3) The naive-vs-corrected comparison gives the paper a concrete methodological payload.

**Concerns at the editorial level.** (1) *Global-change framing is underdeveloped.* The Introduction does not connect to the central global-change debate (net loss vs reorganization; Dornelas, Blowes, McGill, Díaz) or to climate/land-use drivers up front. For GCB, the climate-and-land-use story must lead, not appear only in a correlational results paragraph. (2) *Evidence status.* The self-declared single-chain status (manuscript line 5) and missing credible intervals would, by themselves, prevent a positive editorial decision. (3) *Significance beyond China.* The paper must state explicitly what the result changes for global biodiversity monitoring, not only for China.

**Recommendation:** Major revision. The science can reach GCB; the framing and evidence reporting must be rebuilt.

---

### Reviewer 1 — Methodology

**Overall.** Model choice (spatial multispecies dynamic occupancy via stMsPGOcc; temporal tMsPGOcc breadth extension) is appropriate and current. Execution and reporting are not yet adequate.

**Major issues (what, where, fix).**
1. *Single-chain inference, no diagnostics (Results "Modelled scope"; Methods).* The summary tables report one chain and the manuscript carries no R̂/ESS/WAIC/PPC. **Fix:** rerun ≥4 chains; report max R̂, min bulk/tail ESS, WAIC, and posterior predictive checks for both models; archive fit objects. Until then no estimate is defensible.
2. *Detection–occupancy identifiability under growing effort (core to the paper).* The detection model must demonstrably absorb the temporal growth of effort. **Fix:** report the posterior detection slope on effort and show the occupancy trend is stable under a fixed-effort counterfactual (the manuscript now promises this; it must be delivered with numbers).
3. *Spatial prior scaling.* `phi.unif = c(0.1, 10)` at a 100-km grid may truncate the spatial range. **Fix:** prior-sensitivity analysis; confirm the posterior does not pile against the bound.
4. *Probabilistic Baselga partitioning is non-standard.* Defining expected shared occupancy as the sum of per-species pairwise minima needs justification. **Fix:** derive bounds and decomposition properties in SI; show the partition behaves sensibly under null occupancy fields.
5. *Trend statistics.* Five-period series with OLS is fragile. **Fix:** Theil–Sen + Mann–Kendall (manuscript now states this) and report per-grid uncertainty.

**Minor.** Relaxed closure (5-yr periods) needs the promised 3-yr sensitivity; synchrony = 1.0 is a likely artefact and should be dropped; classification thresholds for expanding/stable/contracting must be stated.

**Scores (0–100):** internal validity 55 (now) → 80 (post-rerun); reproducibility 45 → 85; statistical reporting 40 → 85.

---

### Reviewer 2 — Domain (homogenization & functional diversity)

**Overall.** The biological story is attractive but currently over-claims relative to the evidence shown.

**Major issues.**
1. *Functional change is tiny and uncertainty-free.* Trait volume 1.383→1.365 (−1.3%) and Rao's Q −1.5% are presented without credible intervals. The title-level claim of "functional homogenization" cannot rest on point estimates this small. **Fix:** report posterior CIs and P(decline) for every functional metric; tie the wording strength to the CIs.
2. *Homogenization needs a compositional/spatial signal.* Homogenization is classically defined as rising among-site similarity / falling spatial beta (McKinney & Lockwood; Olden & Rooney). The manuscript leans on temporal beta and alpha-functional means. **Fix:** add the direct test — spatial beta diversity among grids declining through time (the manuscript now promises this; it is essential, not optional).
3. *Literature dialogue is thin for a domain claim.* Must engage Clavel (2010), Devictor (2007, 2012), Le Viol (2012), Jarzyna & Jetz (2017, 2018), and the net-loss-vs-reorganization debate. (The revised draft adds these.)

**Minor.** The winners (waterbirds, open-country, human-associated) vs losers (montane forest specialists) contrast is strong evidence and should be a figure, not just a species list.

**Scores:** originality 75; domain contribution 70 (now) → 85 (with spatial-beta evidence + CIs).

---

### Reviewer 3 — Perspective (monitoring / EBV / citizen science)

**Overall.** The most transferable contribution is methodological: the demonstration that detection correction must be carried into functional and compositional axes. This deserves to be foregrounded for the global-monitoring community (Essential Biodiversity Variables; Jetz et al. 2019; Isaac et al. 2014).

**Major.**
1. *Generalize the message.* State explicitly that EBV-style biodiversity monitoring built on volunteer data will systematically misread homogenization as recovery unless detection and multidimensionality are jointly modelled. This is the GCB-scale "so what".
2. *Effort confound at policy scale.* The strongest contribution is showing where naive and corrected trends disagree — these are the grids where monitoring would mislead. Make this a decision-relevant output (a map + criterion).

**Minor.** Data-availability under provider terms needs a clear, reproducible statement; an EBV-aligned framing would widen the audience.

**Scores:** practical impact 80; cross-disciplinary reach 75.

---

### Reviewer DA — Devil's Advocate

**Strongest counter-argument (≈250 words).**
The paper's central inference — that Chinese birds genuinely expanded — may be an artefact of the very process it claims to correct. Over 2000–2024, birdwatching participation in China grew by orders of magnitude, concentrated in the eastern, lowland, urbanizing grids that the paper reports as the strongest gainers and as positively correlated with cropland, built land and human footprint. A dynamic occupancy model separates detection from occupancy only insofar as the detection covariates capture the true effort process and the closure/colonization structure is correctly specified. If effort growth is imperfectly modelled — plausible, because the detection model is thin (events, duration, missing-duration; source only now being added) and effort and occupancy trends are spatially collinear — then residual detectability is reassigned to colonization, manufacturing exactly the "expansion" the paper reports. The naive-vs-corrected comparison is used as reassurance, but it cuts the other way: only 2–3% of grids flip direction, i.e. corrected and naive trends are nearly collinear, which is what one expects if correction is incomplete rather than sufficient. Meanwhile the functional "homogenization" signal is ~1% with no credible intervals, so the paper's most novel claim may be statistically indistinguishable from zero. The honest null model is: rising effort inflates occupancy of detectable generalists, richness rises, functional space barely moves because no genuinely new functional types are added, and beta diversity looks nested because the same widespread species accumulate everywhere. That null reproduces every headline result without any true biological change.

**Issue list.**
- **CRITICAL** (identifiability): the effort-confound is not yet defeated with numbers. The four promised controls (detection–effort posterior, effort-saturated subset, fixed-effort counterfactual, source term) must be executed and must survive; until then the central claim is unproven. *Decision-blocking per IRON RULE.*
- **CRITICAL** (statistical reality of homogenization): functional declines lack CIs; if they span zero, the title and abstract overstate.
- **MAJOR** (collinearity argument): the small flip rate is presented as strength but is consistent with incomplete correction; address directly.
- **MAJOR** (nestedness ambiguity): late nestedness is read as convergence-by-gain, but could be uneven retention of a common pool; the winners/losers functional identities must adjudicate.
- **MINOR**: "phylogenetic diversity increased" may simply reflect adding any species to a sparse tree; report standardized effect sizes (ses.PD/ses.MPD) against a null.

**Ignored alternative explanations.** Greening/conservation investment (Ouyang 2016; Bryan 2018) could drive genuine gains — a *competing* causal story the paper neither embraces nor excludes. Range shifts tracking warming (climatic debt; Devictor 2012) could masquerade as expansion at the national envelope.

**Missing stakeholders.** National monitoring agencies who might act on a "recovery" headline; the citizen observers whose effort structure is the confound.

**"So what?" test.** Passes *if and only if* the effort confound is defeated and homogenization is shown with uncertainty; otherwise the paper documents an observation artefact.

---

## Phase 2 — Editorial synthesis & decision

### Consensus across reviewers
- **Strong, GCB-relevant idea** (EIC, R2, R3): expansion-without-differentiation is novel and important.
- **Evidence reporting is below bar** (R1, R2, DA): single-chain, no CIs, no executed confound controls.
- **Global-change framing must lead** (EIC, R3): connect to climate/land-use drivers and the net-loss-vs-reorganization debate from the first paragraph.
- **Homogenization needs the spatial-beta test + CIs** (R2, DA): the classic definition is compositional convergence, not an alpha-functional mean drop.

### Disagreement / arbitration
- R3 sees the methodological message (monitoring) as the headline; R2/EIC see the biological homogenization as the headline. **Arbitration:** GCB rewards the *biological* result framed by the *global-change* significance, with the monitoring message as the transferable implication. The revised draft leads with biology, closes with monitoring.

### Devil's Advocate CRITICAL → decision constraint
Per IRON RULE, with two CRITICAL findings open (effort confound not yet defeated with numbers; functional homogenization lacks CIs), **the decision cannot be Accept and cannot be Minor Revision.**

### EDITORIAL DECISION: **Major Revision** (reconsider after substantive revision; not a guarantee of acceptance)

### Revision Roadmap (prioritized; directly actionable)

**Tier 1 — decision-blocking (must be resolved; require reruns/new analysis)**
1. Execute and report the four effort-confound controls with numbers; the expansion signal must survive all four. *(DA-CRITICAL, R1)*
2. Report posterior CIs and P(decline) for all functional metrics; calibrate homogenization wording to the result. *(DA-CRITICAL, R2)*
3. ≥4-chain rerun with R̂/ESS/WAIC/PPC; archive fit objects. *(R1)*
4. Add the direct homogenization test: among-grid spatial beta diversity declining through time. *(R2, DA)*

**Tier 2 — required for GCB framing/quality**
5. Reframe Introduction around global change (drivers + net-loss-vs-reorganization debate + multidimensional facets + China as hotspot), ending in a sharp scientific question + H1–H5. *(EIC, R3)*
6. Expand the Discussion: reconcile with the global debate, defend the confound, separate competing causal stories (effort vs greening/conservation vs climatic debt), and deliver the monitoring/EBV implication. *(EIC, R2, R3, DA)*
7. Strengthen winners–losers into a functional-identity figure; report standardized PD effect sizes. *(R2, DA-minor)*
8. phi prior-sensitivity; 3-yr closure sensitivity; Theil–Sen/Mann–Kendall trends. *(R1)*

**Tier 3 — polish**
9. Expand literature (climate/land-use drivers, homogenization, citizen-science statistics, EBV, China context). *(R2, EIC)* — done in the revised draft with DOI-verified additions.
10. Reproducibility statement (data terms, code, sessionInfo, archive DOI). *(R3)*

> The companion file `MANUSCRIPT_GCB_v6_20260602.md` implements all writing-level items (Tier 2, Tier 3) and writes Tier-1 results as clearly-marked placeholders to be filled from the rerun — no diagnostic value is fabricated.
