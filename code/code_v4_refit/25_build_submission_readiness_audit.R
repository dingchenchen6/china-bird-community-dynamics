#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

root <- Sys.getenv("BIRD_PROJECT_ROOT", getwd())
out_dir <- file.path(root, "editorial_audit_20260721")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exists_any <- function(paths) any(file.exists(file.path(root, paths)))
read_if <- function(path) {
  full <- file.path(root, path)
  if (!file.exists(full)) return(NULL)
  suppressMessages(readr::read_csv(full, show_col_types = FALSE))
}

survey_path <- file.path(root, "data/derived_v3/survey_history_v3.rds")
survey <- if (file.exists(survey_path)) readRDS(survey_path) else NULL
survey_dim <- if (!is.null(survey)) dim(survey$y) else NULL

diag200 <- read_if("results_v3/table_convergence_diagnostics_v3_full_200sp_ar1_spatial_4chain.csv")
diag500 <- read_if("results_v3/table_convergence_diagnostics_v3_full_500sp_ar1_temporal_4chain.csv")
homog <- read_if("results_v3/table_homogenization_trend_v3_full_200sp_ar1_spatial.csv")
effort <- read_if("results_v3/table_effort_confound_controls_v3_full_200sp_ar1_spatial.csv")
manuscript_path <- file.path(root, "results_v4/MANUSCRIPT_GCB_v6_20260602.md")
manuscript <- if (file.exists(manuscript_path)) paste(readLines(manuscript_path, warn = FALSE), collapse = "\n") else ""

checks <- tribble(
  ~domain, ~check, ~status, ~evidence, ~required_action,
  "data", "Visit-level repeated detections available", "PASS", "data/derived_v2 visit_effort and species_visit are present", "Archive these exact inputs with checksums",
  "data", "Canonical survey history is 4D", if (length(survey_dim) == 4L) "PASS" else "FAIL", paste(survey_dim, collapse = " x "), "Rebuild or resync species x site x period x replicate survey history",
  "model", "200-species spatial full fit locally auditable", if (exists_any(c("data/derived_v3/stMsPGOcc_fit_v3_full_200sp_ar1_spatial.rds", "data/derived_v4/stMsPGOcc_fit_v4_full_200sp_ar1_spatial.qs"))) "PASS" else "FAIL", "No full fit object found locally", "Resync all four chain objects and the combined posterior",
  "model", "500-species temporal full fit locally auditable", if (exists_any(c("data/derived_v3/tMsPGOcc_fit_v3_full_500sp_ar1_temporal.rds", "data/derived_v4/tMsPGOcc_fit_v4_full_500sp_ar1_temporal.qs"))) "PASS" else "FAIL", "No full fit object found locally", "Resync all four chain objects and the combined posterior",
  "diagnostics", "200-species monitored parameters meet Rhat < 1.05", if (!is.null(diag200) && max(diag200$rhat, na.rm = TRUE) < 1.05) "PASS" else "FAIL", if (!is.null(diag200)) sprintf("max Rhat %.3f; min ESS %.0f", max(diag200$rhat, na.rm = TRUE), min(diag200$ess, na.rm = TRUE)) else "missing", "Report low-ESS hyperparameters and extend chains if they affect conclusions",
  "diagnostics", "500-species Rhat available", if (!is.null(diag500) && any(is.finite(diag500$rhat))) "PASS" else "FAIL", "All exported Rhat values are NA", "Recover chain identities and recompute split rank-normalized Rhat",
  "diagnostics", "Posterior predictive checks completed", if (exists_any(c("results_v3/table_ppc_bayesian_pvalue.csv", "results_v4/table_ppc_bayesian_pvalue.csv"))) "PASS" else "FAIL", "No PPC result table", "Run detection- and occupancy-level PPC by period and effort stratum",
  "detection", "Data-source detection term fitted", if (!is.null(effort) && any(effort$control == "data_source_detection_term" & is.finite(effort$mean))) "PASS" else "FAIL", "delta WAIC is NA", "Add source/platform composition to the detection model and compare WAIC/trajectory",
  "detection", "Effort-stable subset supports direction", if (!is.null(effort) && any(effort$control == "effort_saturated_subset" & effort$q025 > 0, na.rm = TRUE)) "PASS" else "FAIL", "30-grid subset; positive interval", "Retain as sensitivity only and report its small coverage",
  "homogenization", "Direct among-grid spatial beta trend available", if (!is.null(homog) && any(is.finite(homog$mean_sorensen))) "PASS" else "FAIL", "All five spatial beta estimates are NA", "Recompute from 4D psi posterior before using homogenization in title",
  "phylogeny", "Strict probability-weighted PD reflected in result tables", "RERUN", "Code corrected on 2026-07-21; existing McTavish PD tables predate correction", "Recompute all PD, trends, figures and driver analyses",
  "temporal", "Loreau-de Mazancourt synchrony reflected in result tables", "RERUN", "Denominator corrected on 2026-07-21; existing synchrony values are invalid", "Recompute synchrony and remove old values",
  "sensitivity", "Three-year closure sensitivity completed", if (exists_any(c("results_v3/table_sensitivity_3yr_vs_5yr.csv", "results_v3/table_sensitivity_3yr_vs_5yr_corrected.csv"))) "PASS" else "FAIL", "No completed 3-year output", "Run 3-year and preferably annual-repeat sensitivity",
  "drivers", "Spatially blocked predictive validation", "FAIL", "RF importance is in-sample permutation importance", "Use spatial block cross-validation or present as descriptive association",
  "figures", "Submission figure package with source data", if (dir.exists(file.path(root, "submission_figure_package_20260714_v3"))) "PASS" else "FAIL", "Six-figure package with PNG/PDF/TIFF/PPTX/source data", "Update PD, synchrony and homogenization panels after rerun"
)

write_csv(checks, file.path(out_dir, "submission_readiness_checks.csv"))

claims <- tribble(
  ~claim_id, ~candidate_claim, ~status, ~evidence, ~editorial_wording,
  "C1", "Detection-corrected local richness increased by about one quarter", "SUPPORTED_WITH_CAVEAT", "200- and 500-species summaries agree; detection source/PPC remain incomplete", "Occupancy-derived expected richness increased strongly in both fitted species pools, subject to unresolved observation-process sensitivities.",
  "C2", "Functional diversity contracted", "PARTIAL", "Trait volume and Rao Q point estimates decline; global trend intervals are incomplete", "Functional trait space did not expand with richness and point estimates declined slightly.",
  "C3", "Chinese bird communities homogenized spatially", "NOT_YET_SUPPORTED", "Direct among-grid spatial beta table is all NA", "Recent temporal change became nestedness-dominated, a pattern consistent with convergence but not a direct spatial-homogenization test.",
  "C4", "Phylogenetic diversity increased", "RERUN_REQUIRED", "Existing McTavish PD used an invalid approximation", "Withhold quantitative PD claims until corrected probability-weighted PD is recomputed.",
  "C5", "Generalists drove the expansion", "ASSOCIATIONAL", "Habitat-breadth Spearman signal is positive; phylogenetic regression interval crosses zero", "Expansion was associated with broader habitat use in complete-case analyses.",
  "C6", "Climate and land use drove the changes", "ASSOCIATIONAL", "Variance partitioning and in-sample RF; large spatial fraction and residual", "Community trends covaried with climate and land-use gradients, without causal attribution.",
  "C7", "The 500-species analysis confirms spatial patterns", "UNSUPPORTED", "500-species tMsPGOcc is explicitly non-spatial", "The 500-species model tests breadth of temporal patterns only.",
  "C8", "The analysis is fully reproducible", "NOT_YET_SUPPORTED", "Local 4D survey history and full chain objects are missing", "Code and derived tables are available, but the complete posterior workflow requires resynchronization and rerun."
)
write_csv(claims, file.path(out_dir, "claim_evidence_matrix.csv"))

div500 <- read_if("results_v3/table_diversity_summary_v3_full_500sp_ar1_temporal_extended.csv")
if (!is.null(div500)) {
  metric_status <- div500 |>
    group_by(metric) |>
    summarise(
      n_rows = n(),
      n_finite = sum(is.finite(mean)),
      p1_mean = mean(mean[period == "P1"], na.rm = TRUE),
      p5_mean = mean(mean[period == "P5"], na.rm = TRUE),
      pct_change = 100 * (p5_mean / p1_mean - 1),
      .groups = "drop"
    ) |>
    mutate(
      status = case_when(
        n_finite == 0 ~ "UNUSABLE_ALL_NA",
        metric == "pd_prob_mctavish" ~ "RERUN_PD_ALGORITHM",
        TRUE ~ "AVAILABLE"
      )
    )
  write_csv(metric_status, file.path(out_dir, "metric_status_500sp.csv"))
}

pending_count <- str_count(manuscript, fixed("[PENDING"))
summary_lines <- c(
  "# Automated submission-readiness snapshot",
  "",
  sprintf("Generated: %s", Sys.Date()),
  sprintf("Canonical local survey-history dimensions: %s", paste(survey_dim, collapse = " x ")),
  sprintf("Manuscript PENDING markers: %d", pending_count),
  sprintf("Checks: %d PASS, %d FAIL, %d RERUN", sum(checks$status == "PASS"), sum(checks$status == "FAIL"), sum(checks$status == "RERUN")),
  "",
  "This file is a machine-generated snapshot. See the accompanying editorial review for interpretation."
)
writeLines(summary_lines, file.path(out_dir, "AUTOMATED_AUDIT_SUMMARY.md"))

message("Audit written to: ", out_dir)
