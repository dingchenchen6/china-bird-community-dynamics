#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把占域后验 psi 一致地传播到分类、系统发育、功能、时间动态多样性，
#   不要因二值阈值化或 picante::pd 失去占域校正含义（C4/C12 修复）
#
# Objective / 分析目标:
#   - div_taxonomic       : sum(psi) corrected richness（修复 v2 #12）
#   - div_functional      : 概率加权 CWM/CWSD/trait_volume/Rao Q
#   - div_phylogenetic_prob : 真正概率加权 PD —— 不再二值化 picante::pd
#                            (PD_prob = Σ_e L_e × (1 − ∏(1 − ψ)))
#   - div_feve / div_fdiv : Villéger 2008 概率加权
#   - div_baselga         : 概率加权 Sørensen/Simpson/Nestedness
#   - div_synchrony / div_variance_ratio / div_temporal_turnover :
#     codyn 风格时间动态指标
#   - prepare_trait_matrix: log10 / no-log 分类（修复 v2 双峰 bug）
#   - theil_sen_slope / mann_kendall_test : 优先 mblm/Kendall，
#     fallback 纯 R（修复 v3 O(n²) for 循环 C8）
#
# Input data / 输入数据:
#   psi 向量、性状矩阵、ape::phylo
#
# Main workflow / 主要流程:
#   纯函数定义，无副作用
#
# Key assumptions / 关键假设:
#   - psi 已经是 [0,1] 概率（不需重新归一）
#   - 性状矩阵行序与 psi 名向量对齐
#   - 系统发育树为 ape::phylo 类
#
# Main packages / 主要包:
#   ape, picante（仅 fallback）, mblm, Kendall（可选）
#
# Output directory / 输出路径:
#   不产出文件
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ─────────────────────────────────────────────────────────────────────
# 1. 分类学多样性（修复 v2 #12 bug）
# ─────────────────────────────────────────────────────────────────────

#' 占域校正分类多样性
#' richness = sum(psi)（期望物种数），而非 sum(psi > eps)
#' shannon / inv_simpson 用归一化 psi 计算
div_taxonomic <- function(psi, eps = 1e-12) {
  psi <- pmax(as.numeric(psi), 0)
  if (sum(psi) < eps) {
    return(list(richness = 0, shannon = 0, inv_simpson = 1))
  }
  richness <- sum(psi)
  p <- psi / sum(psi)
  shannon <- -sum(p[p > eps] * log(p[p > eps]))
  inv_simpson <- 1 / sum(p^2)
  list(richness = richness, shannon = shannon, inv_simpson = inv_simpson)
}

# ─────────────────────────────────────────────────────────────────────
# 2. 功能多样性（概率加权）
# ─────────────────────────────────────────────────────────────────────

#' 概率加权 CWM / CWSD / trait_volume / trait_dispersion / Rao Q
#' v4 改进：eps 默认 PSI_EPS_DEFAULT（0.05），可参数化
div_functional <- function(psi, trait_mat, eps = NULL) {
  if (is.null(eps)) eps <- if (exists("PSI_EPS_DEFAULT")) PSI_EPS_DEFAULT else 1e-12
  ok <- psi > eps
  n_trait <- ncol(trait_mat)

  if (!any(ok) || is.null(trait_mat) || nrow(trait_mat) != length(psi)) {
    return(list(
      cwm = setNames(rep(NA_real_, n_trait), colnames(trait_mat)),
      cwsd = setNames(rep(NA_real_, n_trait), colnames(trait_mat)),
      trait_volume = NA_real_,
      trait_dispersion = NA_real_,
      rao_q = NA_real_
    ))
  }

  psi_ok  <- psi[ok]
  trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)

  cwm <- colSums(w * trait_ok)
  diff_mat <- sweep(trait_ok, 2, cwm, "-")
  cwsd <- sqrt(colSums(w * diff_mat^2))

  cov_mat <- crossprod(diff_mat * sqrt(w))
  trait_volume <- sqrt(sum(diag(cov_mat)))
  trait_dispersion <- sum(w * sqrt(rowSums(diff_mat^2)))

  dist_mat <- as.matrix(stats::dist(trait_ok))
  rao_q <- sum(outer(w, w) * dist_mat)

  list(cwm = cwm, cwsd = cwsd,
       trait_volume = trait_volume,
       trait_dispersion = trait_dispersion,
       rao_q = rao_q)
}

# ─────────────────────────────────────────────────────────────────────
# 3. 系统发育多样性 —— v4 真正概率加权 PD（C4 关键修复）
# ─────────────────────────────────────────────────────────────────────

#' 真正概率加权 Faith's PD
#' 算法：对系统发育树每条边 e（枝长 L_e）：
#'   该边的"存在概率" = 1 − ∏_{sp ∈ desc(e)}(1 − ψ_sp)
#'   PD_prob = Σ_e L_e × (存在概率)
#' 这是 manuscript 描述的公式；v3 之前用 picante::pd 二值化是不一致的。
#'
#' @param psi 物种占有率向量（命名或 tip_order 对齐）
#' @param phylo ape::phylo
#' @param tip_order character：物种顺序（与 psi 对齐）；NULL 时用 names(psi)
#' @param eps 数值精度（用于判断"存在")
#' @return list: pd_prob, mpd_prob
div_phylogenetic_prob <- function(psi, phylo, tip_order = NULL, eps = 1e-6) {
  if (is.null(phylo)) return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  if (!inherits(phylo, "phylo")) return(list(pd_prob = NA_real_, mpd_prob = NA_real_))

  psi <- pmax(as.numeric(psi), 0)
  if (!is.null(tip_order)) names(psi) <- tip_order
  if (is.null(names(psi))) return(list(pd_prob = NA_real_, mpd_prob = NA_real_))

  matched <- intersect(names(psi), phylo$tip.label)
  if (length(matched) < 2) return(list(pd_prob = NA_real_, mpd_prob = NA_real_))

  psi_m <- psi[matched]
  if (sum(psi_m) < eps) return(list(pd_prob = NA_real_, mpd_prob = NA_real_))

  # 剪枝
  pruned <- ape::drop.tip(phylo, setdiff(phylo$tip.label, matched))
  tips <- pruned$tip.label
  psi_t <- psi_m[tips]

  # ── 对每条边：找其下游所有 tip → 计算"边存在概率" ───────────────────
  n_tip <- length(tips)
  edges <- pruned$edge          # 2 列：parent, child
  edge_lengths <- pruned$edge.length

  # 对每条边的 child 节点，找其下游所有 tip 的 indices
  # 用 ape::prop.part 或自己 BFS；这里用 phangorn 风格 Descendants
  descendants <- .descendants_per_edge(pruned)

  # 每条边的"被任意后代物种占据"的概率 = 1 − ∏_{sp ∈ desc}(1 − ψ_sp)
  edge_pres <- sapply(descendants, function(tip_idx) {
    if (length(tip_idx) == 0) return(0)
    psi_d <- psi_t[tip_idx]
    1 - prod(1 - psi_d)
  })

  pd_prob <- sum(edge_lengths * edge_pres)

  # ── MPD（概率加权平均系统发育距离） ────────────────────────────────
  cophen <- ape::cophenetic.phylo(pruned)
  w <- psi_t / sum(psi_t)
  w_mat <- outer(w, w); diag(w_mat) <- 0
  mpd_prob <- if (sum(w_mat) > 0) sum(cophen * w_mat) / sum(w_mat) else NA_real_

  list(pd_prob = pd_prob, mpd_prob = mpd_prob)
}

# 内部：每条边的后代 tip indices
.descendants_per_edge <- function(phylo) {
  n_tip <- length(phylo$tip.label)
  n_node <- phylo$Nnode
  edges <- phylo$edge

  # 对每个 internal node 缓存其后代 tip indices
  # tip nodes 编号 1..n_tip，internal node 编号 (n_tip+1)..(n_tip+n_node)
  child_list <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(edges))) {
    p <- edges[i, 1]; c <- edges[i, 2]
    child_list[[p]] <- c(child_list[[p]], c)
  }

  # 递归：node 的后代 tips
  desc_cache <- vector("list", n_tip + n_node)
  get_desc <- function(node) {
    if (!is.null(desc_cache[[node]])) return(desc_cache[[node]])
    if (node <= n_tip) {
      desc_cache[[node]] <<- node
      return(node)
    }
    children <- child_list[[node]]
    out <- unlist(lapply(children, get_desc))
    desc_cache[[node]] <<- out
    out
  }

  # 对每条边的 child 节点取后代
  lapply(seq_len(nrow(edges)), function(i) get_desc(edges[i, 2]))
}

# 兼容旧名：div_phylogenetic 默认走概率加权版本
div_phylogenetic <- div_phylogenetic_prob

# ─────────────────────────────────────────────────────────────────────
# 4. 功能均匀度 / 分散度 FEve / FDiv（Villéger et al. 2008）
# ─────────────────────────────────────────────────────────────────────

div_feve <- function(psi, trait_mat, eps = NULL) {
  if (is.null(eps)) eps <- if (exists("PSI_EPS_DEFAULT")) PSI_EPS_DEFAULT else 1e-12
  ok <- psi > eps
  S <- sum(ok)
  if (S < 3 || is.null(trait_mat) || nrow(trait_mat) != length(psi)) return(NA_real_)

  psi_ok <- psi[ok]; trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)
  D <- as.matrix(stats::dist(trait_ok))

  # Prim MST
  n <- S
  in_tree <- rep(FALSE, n); min_dist <- rep(Inf, n); min_from <- rep(NA_integer_, n)
  min_dist[1] <- 0
  edge_i <- integer(n - 1); edge_j <- integer(n - 1)

  for (k in seq_len(n)) {
    u <- which.min(ifelse(in_tree, Inf, min_dist))
    if (is.na(u) || min_dist[u] == Inf) break
    in_tree[u] <- TRUE
    if (k > 1) { edge_i[k - 1] <- min_from[u]; edge_j[k - 1] <- u }
    for (v in seq_len(n)) {
      if (!in_tree[v] && D[u, v] < min_dist[v]) {
        min_dist[v] <- D[u, v]; min_from[v] <- u
      }
    }
  }
  n_edges <- sum(edge_i > 0)
  if (n_edges < 2) return(NA_real_)

  EW <- numeric(n_edges)
  for (k in seq_len(n_edges)) {
    EW[k] <- (w[edge_i[k]] + w[edge_j[k]]) * D[edge_i[k], edge_j[k]] / 2
  }
  sum_EW <- sum(EW)
  if (sum_EW < eps) return(NA_real_)
  PEW <- EW / sum_EW
  PEW_ref <- 1 / (S - 1)
  FEve <- 1 - sum(abs(PEW - PEW_ref)) / (2 * (1 - PEW_ref))
  max(0, min(1, FEve))
}

div_fdiv <- function(psi, trait_mat, eps = NULL) {
  if (is.null(eps)) eps <- if (exists("PSI_EPS_DEFAULT")) PSI_EPS_DEFAULT else 1e-12
  ok <- psi > eps
  if (sum(ok) < 3 || is.null(trait_mat) || nrow(trait_mat) != length(psi)) return(NA_real_)
  psi_ok <- psi[ok]; trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)
  G <- colSums(w * trait_ok)
  dG <- sqrt(rowSums(sweep(trait_ok, 2, G, "-")^2))
  dbar <- sum(w * dG); dmax <- max(dG)
  if (dmax < eps) return(NA_real_)
  num <- sum(w * abs(dG - dbar)) + dbar
  den <- sum(w * dG) + dbar
  if (den < eps) return(NA_real_)
  fdiv <- num / den
  max(0, min(1, fdiv))
}

# ─────────────────────────────────────────────────────────────────────
# 5. Baselga 时间 β 分解（概率加权）
# ─────────────────────────────────────────────────────────────────────

div_baselga <- function(psi_t1, psi_t2) {
  psi_t1 <- pmax(as.numeric(psi_t1), 0)
  psi_t2 <- pmax(as.numeric(psi_t2), 0)

  a <- sum(pmin(psi_t1, psi_t2))
  b <- max(0, sum(psi_t1) - a)
  c <- max(0, sum(psi_t2) - a)

  denom_sor <- 2 * a + b + c
  denom_sim <- a + min(b, c)

  beta_sor <- if (denom_sor > 0) (b + c) / denom_sor else NA_real_
  beta_sim <- if (denom_sim > 0) min(b, c) / denom_sim else NA_real_
  beta_sne <- if (!is.na(beta_sor) && !is.na(beta_sim)) beta_sor - beta_sim else NA_real_

  prop_turnover <- if (!is.na(beta_sor) && beta_sor > 1e-10) {
    min(1, max(0, beta_sim / beta_sor))
  } else NA_real_

  list(beta_sor = beta_sor, beta_sim = beta_sim, beta_sne = beta_sne,
       a = a, b = b, c = c, prop_turnover = prop_turnover)
}

# ─────────────────────────────────────────────────────────────────────
# 6. 时间动态：synchrony / variance_ratio / turnover
# ─────────────────────────────────────────────────────────────────────

div_synchrony <- function(comm_matrix, eps = 1e-12) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2) return(NA_real_)
  sp_var <- apply(comm_matrix, 1, var, na.rm = TRUE)
  keep <- sp_var > eps
  if (sum(keep) < 2) return(NA_real_)
  comm_matrix <- comm_matrix[keep, , drop = FALSE]
  total_ts  <- colSums(comm_matrix)
  var_total <- var(total_ts, na.rm = TRUE)
  sum_sd_sp <- sum(apply(comm_matrix, 1, sd, na.rm = TRUE))
  if (sum_sd_sp < eps) return(NA_real_)
  synchrony <- var_total / sum_sd_sp^2
  min(1, max(0, synchrony))
}

div_variance_ratio <- function(comm_matrix, eps = 1e-12) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2) return(NA_real_)
  sp_var <- apply(comm_matrix, 1, var, na.rm = TRUE)
  keep <- sp_var > eps
  if (sum(keep) < 2) return(NA_real_)
  comm_matrix <- comm_matrix[keep, , drop = FALSE]
  total_ts  <- colSums(comm_matrix)
  var_total <- var(total_ts, na.rm = TRUE)
  sum_var_sp <- sum(apply(comm_matrix, 1, var, na.rm = TRUE))
  if (sum_var_sp < eps) return(NA_real_)
  var_total / sum_var_sp
}

div_temporal_turnover <- function(comm_matrix, threshold = 0.1) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2)
    return(list(turnover_total = NA_real_, turnover_gain = NA_real_, turnover_loss = NA_real_))
  n_periods <- ncol(comm_matrix)
  total_gain <- 0; total_loss <- 0; n_pairs <- 0
  for (t in seq_len(n_periods - 1)) {
    present_t  <- comm_matrix[, t]     > threshold
    present_t1 <- comm_matrix[, t + 1] > threshold
    gain <- sum(!present_t & present_t1)
    loss <- sum(present_t & !present_t1)
    pool <- sum(present_t | present_t1)
    if (pool > 0) {
      total_gain <- total_gain + gain / pool
      total_loss <- total_loss + loss / pool
      n_pairs <- n_pairs + 1
    }
  }
  if (n_pairs == 0)
    return(list(turnover_total = NA_real_, turnover_gain = NA_real_, turnover_loss = NA_real_))
  list(turnover_total = (total_gain + total_loss) / n_pairs,
       turnover_gain  = total_gain / n_pairs,
       turnover_loss  = total_loss / n_pairs)
}

# ─────────────────────────────────────────────────────────────────────
# 7. 性状矩阵预处理（v3 修复 #13 沿用）
# ─────────────────────────────────────────────────────────────────────

prepare_trait_matrix <- function(trait_aligned, trait_vars,
                                  log_vars = TRAIT_VARS_LOG10) {
  no_log_vars <- setdiff(trait_vars, log_vars)

  log_part <- if (length(intersect(log_vars, trait_vars)) > 0) {
    trait_aligned |>
      select(any_of(intersect(log_vars, trait_vars))) |>
      mutate(across(where(is.numeric),
                    ~ if_else(.x > 0, log10(.x), NA_real_)))
  } else {
    tibble::tibble(.rows = nrow(trait_aligned))
  }

  no_log_part <- if (length(intersect(no_log_vars, trait_vars)) > 0) {
    trait_aligned |> select(any_of(intersect(no_log_vars, trait_vars)))
  } else {
    tibble::tibble(.rows = nrow(trait_aligned))
  }

  trait_mat <- bind_cols(log_part, no_log_part) |>
    select(all_of(trait_vars)) |>
    as.matrix()
  rownames(trait_mat) <- trait_aligned$species

  for (j in seq_len(ncol(trait_mat))) {
    miss <- is.na(trait_mat[, j])
    if (any(miss)) trait_mat[miss, j] <- median(trait_mat[, j], na.rm = TRUE)
  }

  trait_mat <- scale(trait_mat)
  attr(trait_mat, "scaled:center") <- NULL
  attr(trait_mat, "scaled:scale")  <- NULL
  trait_mat
}

# ─────────────────────────────────────────────────────────────────────
# 8. 稳健趋势估计（C8 修复：优先 mblm / Kendall，fallback 纯 R）
# ─────────────────────────────────────────────────────────────────────

#' Theil-Sen 斜率（优先 mblm；fallback 纯 R）
theil_sen_slope <- function(y, x = seq_along(y)) {
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < 3) return(NA_real_)

  if (requireNamespace("mblm", quietly = TRUE)) {
    df <- data.frame(y = y, x = x)
    fit <- tryCatch(mblm::mblm(y ~ x, dataframe = df, repeated = FALSE),
                     error = function(e) NULL)
    if (!is.null(fit)) return(unname(coef(fit)[2]))
  }

  # fallback: 纯 R O(n²)，n=5 时仍快
  slopes <- numeric(n * (n - 1) / 2)
  k <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      slopes[k] <- (y[j] - y[i]) / (x[j] - x[i])
    }
  }
  median(slopes[1:k])
}

#' Mann-Kendall 趋势检验（优先 Kendall::MannKendall；fallback 纯 R）
mann_kendall_test <- function(y) {
  ok <- !is.na(y)
  y <- y[ok]; n <- length(y)
  if (n < 4) return(list(S = NA_real_, tau = NA_real_, p_value = NA_real_))

  if (requireNamespace("Kendall", quietly = TRUE)) {
    mk <- tryCatch(Kendall::MannKendall(y), error = function(e) NULL)
    if (!is.null(mk)) {
      return(list(
        S = unname(mk$S),
        tau = unname(mk$tau),
        p_value = unname(mk$sl)
      ))
    }
  }

  # fallback：纯 R
  S <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) S <- S + sign(y[j] - y[i])
  }
  tau <- S / (n * (n - 1) / 2)
  tie_groups <- table(y)
  tie_correction <- sum(tie_groups * (tie_groups - 1) * (2 * tie_groups + 5))
  var_S <- (n * (n - 1) * (2 * n + 5) - tie_correction) / 18
  Z <- if (S > 0) (S - 1) / sqrt(max(var_S, 1)) else if (S < 0) (S + 1) / sqrt(max(var_S, 1)) else 0
  list(S = S, tau = tau, p_value = 2 * pnorm(-abs(Z)))
}

# ─────────────────────────────────────────────────────────────────────
# 9. 后验汇总辅助
# ─────────────────────────────────────────────────────────────────────

#' 对 draws × ... 数组做按位置后验汇总
summarise_post <- function(arr, probs = c(0.025, 0.5, 0.975)) {
  if (is.null(arr)) return(tibble::tibble())
  if (is.vector(arr)) arr <- matrix(arr, ncol = 1)
  out <- apply(arr, length(dim(arr)), function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd   = sd(x, na.rm = TRUE),
      q025 = quantile(x, probs[1], na.rm = TRUE),
      q50  = quantile(x, probs[2], na.rm = TRUE),
      q975 = quantile(x, probs[3], na.rm = TRUE))
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  tibble::as_tibble(t(out))
}

message("[utils_diversity_v4] loaded — probabilistic PD, mblm/Kendall, eps configurable")
