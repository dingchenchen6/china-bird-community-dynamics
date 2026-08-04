#!/usr/bin/env Rscript
# ============================================================
# 32_conservation_prioritization.R
#
# Scientific question / 科学问题:
#   如果保护规划以物种数为目标，会不会系统性错过"阻止功能同质化"最关键的区域？
#   本脚本用系统保护规划（最小集覆盖 / 最大效用），对比三套优先区方案：
#     (A) 物种数导向 richness-led
#     (B) 功能独特性导向 functional-distinctiveness-led
#     (C) 同质化风险导向 homogenization-risk-led
#   并量化三者的空间错配，以及现有保护区对各方案的满足程度。
#   Does richness-led planning systematically miss the areas that matter
#   for preventing functional homogenization?
#
# Objective / 分析目标:
#   1) 构建三类保护特征层（物种、功能独特性加权、同质化风险）
#   2) 以 prioritizr 求解 30% 陆域预算下的最优保护网络（对标昆明-蒙特利尔 30x30）
#   3) 量化方案间空间重叠（Jaccard / Cohen's kappa）——错配即政策含义
#   4) 评估现有保护区对三套方案的达成率，识别扩建优先级
#
# Input / 输入:
#   derived : psi_samples_thinned_*  (4D)   物种占域后验
#   results : table_pa_grid_coverage_*      现有保护区覆盖率（脚本31产出）
#   results : table_trend_summary_*         功能体积趋势（同质化风险）
#   results : table_species_trend_traits_*  性状（功能独特性）
# Output / 输出:
#   results: table_priority_solutions_<run_label>.csv    三方案逐网格选中情况
#   results: table_priority_overlap_<run_label>.csv      方案间重叠度
#   results: table_priority_pa_shortfall_<run_label>.csv 现有保护区达成率与缺口
#
# Key assumptions / 关键假设:
#   - 成本层用等面积（每网格成本相同）；若有土地机会成本数据应替换
#   - 功能独特性用性状空间中的平均距离（越独特权重越高）
#   - 30% 预算对标 Kunming-Montreal GBF Target 3；另跑 17% 作敏感性
# Packages: prioritizr, dplyr, readr (求解器: highs 或 Rsymphony)
# References: Margules & Pressey 2000 Nature; Rodrigues et al. 2004 Nature;
#             Pollock et al. 2017 Nature; Brum et al. 2017 PNAS;
#             Jung et al. 2021 Nat Ecol Evol; Hanson et al. 2024 Conserv Biol
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

ensure_v3_dirs()
is_pilot  <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
BUDGETS <- c(0.30, 0.17)     # GBF Target 3 与旧 Aichi 目标
log_time("32", "Systematic conservation prioritization")

# ── 1. 载入占域后验与辅助层 ─────────────────────────────────────────
psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
if (is.null(psi_obj) || length(dim(psi_obj$psi_samples_thinned)) < 4)
  stop("[32] need 4D psi_samples_thinned; re-thin via 05.")
psi <- apply(psi_obj$psi_samples_thinned, c(2, 3, 4), mean)   # sp × site × period
sp <- psi_obj$species; sites <- psi_obj$sites
np <- dim(psi)[3]; ns <- length(sp); ng <- length(sites)
feat <- psi[, , np]                                           # 用最新期占域作为特征层

cov_path <- paste0(v3_file("results", paste0("table_pa_grid_coverage_", run_label)), ".csv")
pa_frac <- if (file.exists(cov_path)) {
  cv <- read_csv(cov_path, show_col_types = FALSE)
  x <- cv$pa_frac_all[match(sites, cv$grid_cell)]; x[is.na(x)] <- 0; x
} else { message("[32] 无保护区覆盖表，达成率分析将跳过"); rep(NA_real_, ng) }

# ── 2. 特征权重：物种数导向 vs 功能独特性导向 ───────────────────────
# 功能独特性 = 该物种在性状空间中到其它物种的平均距离（越大越独特）
tr_path <- paste0(v3_file("results", paste0("table_species_trend_traits_", run_label)), "_extended.csv")
w_rich <- rep(1, ns)                                          # (A) 等权 = 物种数导向
w_func <- rep(1, ns)
if (file.exists(tr_path)) {
  tr <- read_csv(tr_path, show_col_types = FALSE)
  # 优先使用 AVONET 形态性状空间（与稿件 §2.6 的功能多样性口径一致）；
  # 缺失时退回生态属性评分，并在日志中标明，避免"功能空间"名实不符。
  # Prefer the AVONET morphological space used elsewhere in the paper.
  morph <- intersect(c("Beak.Length_Culmen", "Beak.Width", "Beak.Depth",
                       "Wing.Length", "Tarsus.Length", "Tail.Length", "Mass",
                       "beak_length_culmen", "wing_length", "tarsus_length", "mass"),
                     names(tr))
  tcols <- if (length(morph) >= 3) morph else
    intersect(c("habitat_breadth", "diet_specialization", "migration_score"), names(tr))
  message("[32] 功能独特性性状空间：", if (length(morph) >= 3)
          paste0("AVONET 形态（", length(morph), " 项）") else "生态属性评分（形态性状缺失）")
  if (length(tcols) >= 2) {
    M <- tr[match(sp, tr$species), tcols, drop = FALSE] |> as.matrix()
    M <- scale(M)
    M[is.na(M)] <- 0
    d <- as.matrix(dist(M))
    uniq <- rowMeans(d)                                       # 平均性状距离
    w_func <- as.numeric(uniq / mean(uniq))                   # (B) 功能独特性权重
    message(sprintf("[32] 功能独特性权重: 范围 %.2f–%.2f", min(w_func), max(w_func)))
  }
}

# (C) 同质化风险导向：功能体积下降越快的网格，保护价值越高
trend_path <- paste0(v3_file("results", paste0("table_trend_summary_", run_label)), "_extended.csv")
risk <- rep(1, ng)
if (file.exists(trend_path)) {
  tt <- read_csv(trend_path, show_col_types = FALSE) |>
    filter(metric %in% c("fric_prob", "trait_volume")) |>
    group_by(grid_cell) |> summarise(slope = mean(mean, na.rm = TRUE), .groups = "drop")
  s <- tt$slope[match(sites, tt$grid_cell)]
  s[is.na(s)] <- median(s, na.rm = TRUE)
  risk <- as.numeric(scale(-s)); risk <- risk - min(risk) + 0.1   # 下降越快风险越高
  message(sprintf("[32] 同质化风险层: 范围 %.2f–%.2f", min(risk), max(risk)))
}

# ── 3. 求解三套方案 ─────────────────────────────────────────────────
if (!requireNamespace("prioritizr", quietly = TRUE)) {
  message("[32] prioritizr 未安装；install.packages(c('prioritizr','highs'))")
  quit(save = "no", status = 0)
}
library(prioritizr)

solve_plan <- function(feature_mat, budget_frac, label) {
  pu <- data.frame(id = seq_len(ng), cost = 1, locked_in = FALSE)
  n_sel <- floor(budget_frac * ng)
  p <- problem(pu, features = data.frame(id = seq_len(nrow(feature_mat)),
                                         name = rownames(feature_mat)),
               rij_matrix = feature_mat, cost_column = "cost") |>
    add_max_utility_objective(budget = n_sel) |>
    add_binary_decisions() |>
    add_default_solver(gap = 0.05, verbose = FALSE)
  s <- solve(p)
  data.frame(grid_cell = sites, selected = s$solution_1,
             plan = label, budget = budget_frac, stringsAsFactors = FALSE)
}

# ⚠️ 方案 B 的构造已修正。原实现仅用物种层权重 w_func 加权：由于站点效用是
#    200 个物种的加总，而 w_func 的变异系数仅约 0.15，加总后站点排序几乎不变，
#    导致方案 B 与 A 的 Jaccard 高达 0.96 —— 那是方法学假象，不是生态结论。
#    方案 C 之所以不同，是因为它用的是站点层权重(risk)，直接改变各网格价值。
#    公平比较要求 B 也在站点层设定目标：这里用"站点功能独特性"——
#    网格内物种的平均性状独特性（占域概率加权），再与物种层权重相乘。
# Plan B rebuilt: species-level weights alone barely change site ranking
# (CV of w_func ~= 0.15), which made B nearly identical to A by construction.
# Plan C differed only because it used site-level weights. B now also carries a
# site-level functional target so the three plans are structurally comparable.
site_func <- as.numeric((w_func %*% feat) / pmax(colSums(feat), 1e-9))   # 站点平均功能独特性
site_func <- as.numeric(scale(site_func)); site_func <- site_func - min(site_func) + 0.1
message(sprintf("[32] 站点功能独特性层: 范围 %.2f–%.2f（CV %.3f）",
                min(site_func), max(site_func), sd(site_func) / mean(site_func)))
message(sprintf("[32] 物种层功能权重 w_func 的 CV = %.3f（若过小则单靠物种权重无法改变选址）",
                sd(w_func) / mean(w_func)))

feat_rich <- feat * w_rich; rownames(feat_rich) <- sp
feat_func <- sweep(feat * w_func, 2, site_func, "*"); rownames(feat_func) <- sp
feat_risk <- sweep(feat, 2, risk, "*"); rownames(feat_risk) <- sp

sols <- do.call(rbind, lapply(BUDGETS, function(b) rbind(
  solve_plan(feat_rich, b, "A_richness_led"),
  solve_plan(feat_func, b, "B_functional_distinctiveness"),
  solve_plan(feat_risk, b, "C_homogenization_risk"))))
write_csv(sols, paste0(v3_file("results", paste0("table_priority_solutions_", run_label)), ".csv"))

# ── 4. 方案间空间错配：错配即政策含义 ────────────────────────────────
ov <- do.call(rbind, lapply(BUDGETS, function(b) {
  w <- filter(sols, budget == b)
  m <- tapply(w$selected, list(w$grid_cell, w$plan), max)
  cmb <- combn(colnames(m), 2, simplify = FALSE)
  do.call(rbind, lapply(cmb, function(p) {
    a <- m[, p[1]] == 1; c2 <- m[, p[2]] == 1
    data.frame(budget = b, plan_1 = p[1], plan_2 = p[2],
               jaccard = sum(a & c2) / sum(a | c2),
               n_only_1 = sum(a & !c2), n_only_2 = sum(!a & c2),
               stringsAsFactors = FALSE)
  }))
}))
print(ov); write_csv(ov, paste0(v3_file("results", paste0("table_priority_overlap_", run_label)), ".csv"))
message("[32] >>> Jaccard 越低 = 物种数导向与功能导向选出的地方越不同，")
message("[32] >>> 即：只按物种数规划会系统性错过阻止同质化的关键区域。")

# ── 5. 现有保护区达成率与扩建缺口 ───────────────────────────────────
if (!all(is.na(pa_frac))) {
  short <- do.call(rbind, lapply(unique(sols$plan), function(pl) {
    w <- filter(sols, plan == pl, budget == 0.30)
    sel <- w$selected == 1
    data.frame(plan = pl,
               n_priority_grids = sum(sel),
               mean_existing_pa_cover = mean(pa_frac[sel]),
               n_priority_wellprotected = sum(pa_frac[sel] >= 0.30),
               pct_priority_wellprotected = 100 * mean(pa_frac[sel] >= 0.30),
               n_priority_unprotected = sum(pa_frac[sel] < 0.05),
               stringsAsFactors = FALSE)
  }))
  print(short)
  write_csv(short, paste0(v3_file("results", paste0("table_priority_pa_shortfall_", run_label)), ".csv"))
  message("[32] >>> pct_priority_wellprotected 低 = 现有保护区网络与优先区错配，")
  message("[32] >>> n_priority_unprotected 即为扩建/新建的首选目标网格。")
}

log_time("32", "DONE")
