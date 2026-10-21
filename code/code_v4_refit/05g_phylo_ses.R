#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题: 系统发育多样性的变化是否超出物种数变化本身所能解释的范围？
#   Does phylogenetic diversity change beyond what richness change alone implies?
# Objective / 分析目标: 用修正后的估计量（枝长按 1 − ∏(1 − ψ) 加权）计算概率加权 PD 与 MPD，
#   并以"打乱物种在系统树上的位置"的零模型计算 ses.PD、ses.MPD。（BLOCKERS Gate 3；SI S5）
# Input / 输入: psi_samples_thinned_<RUN_LABEL>.rds（抽样 × 物种 × 格点 × 时期）；phylogeny_matched_v3.rds
# Output / 输出: results_v3/table_phylo_ses_<RUN_LABEL>.csv（格点 × 时期）、table_phylo_ses_trend_<RUN_LABEL>.csv
# Key assumptions / 关键假设: 用后验均值 ψ 计算（可用 PHYLO_DRAWS 改为逐抽样）；零模型保持各格点 ψ 向量
#   不变、只交换物种在树上的标签（taxa.labels），默认 199 次。
# Packages / 主要包: ape
# ============================================================
suppressPackageStartupMessages(library(ape))
CODE <- Sys.getenv("V3_CODE_DIR", "code_v4_refit")
source(file.path(CODE, "00_config.R")); source(file.path(CODE, "utils_paths.R")); source(file.path(CODE, "utils_diversity_extended.R"))
N_NULL <- as.integer(Sys.getenv("PHYLO_N_NULL", "199")); set.seed(20261020)
ps <- readRDS(v3_file("derived", paste0("psi_samples_thinned_", RUN_LABEL), "rds"))
arr <- ps$psi_samples_thinned; stopifnot(length(dim(arr)) == 4)
tree <- readRDS(v3_file("derived", "phylogeny_matched_v3", "rds"))
if (!inherits(tree, "phylo")) stop("[05g] phylogeny_matched_v3.rds 不是 phylo 对象")
sp <- ps$species; keep <- sp %in% tree$tip.label
message(sprintf("[05g] 系统树匹配 %d/%d 种", sum(keep), length(sp)))
psi_mean <- apply(arr[, keep, , , drop = FALSE], c(2, 3, 4), mean)              # 物种 × 格点 × 时期
spk <- sp[keep]; tree <- keep.tip(tree, spk)

pd_mpd <- function(p, labels) {
  names(p) <- labels
  r <- div_phylogenetic_generic(p[tree$tip.label], tree)
  c(pd = r$pd_prob, mpd = r$mpd_prob)
}
rows <- list()
for (t in seq_len(dim(psi_mean)[3])) for (i in seq_len(dim(psi_mean)[2])) {
  p <- psi_mean[, i, t]; if (sum(p > 1e-6) < 3) next
  obs <- pd_mpd(p, spk)
  nul <- replicate(N_NULL, pd_mpd(p, sample(spk)))                              # 交换标签 / taxa-label shuffle
  ses <- (obs - rowMeans(nul)) / apply(nul, 1, sd)
  rows[[length(rows) + 1]] <- data.frame(grid_cell = ps$sites[i], period = paste0("P", t),
    pd_prob = obs["pd"], mpd_prob = obs["mpd"], ses_pd = ses["pd"], ses_mpd = ses["mpd"])
}
tab <- do.call(rbind, rows)
write.csv(tab, v3_file("results", paste0("table_phylo_ses_", RUN_LABEL)), row.names = FALSE)
tr <- aggregate(cbind(pd_prob, mpd_prob, ses_pd, ses_mpd) ~ period, tab, mean)
write.csv(tr, v3_file("results", paste0("table_phylo_ses_trend_", RUN_LABEL)), row.names = FALSE)
print(tr)
