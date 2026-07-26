## utils_diversity.R
## 给定一个 species × 1 的 psi 向量（一个 grid × period 的占域概率），
## 返回该群落的多样性指标。所有函数都是矢量友好、对 NA 安全的。
##
## 设计目标：让 stage-5 后处理在每一个 MCMC draw 上调用这些函数，把
## 后验不确定性自然传播到群落多样性指标层。

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
})

## ---------- 1. 分类 / 计数多样性 ------------------------------------------

div_taxonomic <- function(psi, eps = 1e-12) {
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps
  ps  <- psi[ok]
  if (length(ps) == 0) {
    return(c(corrected_richness = 0,
             shannon            = NA_real_,
             inverse_simpson    = NA_real_,
             effective_species  = NA_real_))
  }
  S   <- sum(ps)
  pn  <- ps / S
  H   <- -sum(pn * log(pn))
  iS  <- 1 / sum(pn^2)
  c(corrected_richness = S,
    shannon            = H,
    inverse_simpson    = iS,
    effective_species  = exp(H))
}

## ---------- 2. CWM / CWSD / 性状空间 --------------------------------------

div_functional <- function(psi, trait_mat, eps = 1e-12) {
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps
  if (!any(ok) || is.null(trait_mat) || nrow(trait_mat) != length(psi)) {
    return(list(cwm = setNames(rep(NA_real_, ncol(trait_mat)),
                                paste0("cwm_", colnames(trait_mat))),
                cwsd = setNames(rep(NA_real_, ncol(trait_mat)),
                                 paste0("cwsd_", colnames(trait_mat))),
                trait_volume = NA_real_,
                trait_dispersion = NA_real_))
  }
  w <- psi[ok]; w <- w / sum(w)
  T <- trait_mat[ok, , drop = FALSE]
  cwm <- colSums(T * w)
  cwsd <- sqrt(colSums(t(t(T) - cwm)^2 * w))
  # trait volume: |Cov| 行列式（稀疏时返回 NA）
  vol <- if (sum(ok) >= max(2, ncol(T) + 1)) {
    cv <- stats::cov.wt(T, wt = w)$cov
    eg <- eigen(cv, symmetric = TRUE, only.values = TRUE)$values
    if (any(eg <= 1e-12)) NA_real_ else exp(sum(log(eg)))
  } else NA_real_
  # trait dispersion: 加权到质心的距离均值
  disp <- sum(sqrt(rowSums((t(t(T) - cwm))^2)) * w)
  list(
    cwm  = setNames(cwm,  paste0("cwm_",  colnames(T))),
    cwsd = setNames(cwsd, paste0("cwsd_", colnames(T))),
    trait_volume     = vol,
    trait_dispersion = disp
  )
}

## ---------- 3. Rao's Q（基于性状距离矩阵） --------------------------------

div_rao_q <- function(psi, dist_mat, eps = 1e-12) {
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps
  if (sum(ok) < 2) return(NA_real_)
  w <- psi[ok]; w <- w / sum(w)
  D <- as.matrix(dist_mat)[ok, ok, drop = FALSE]
  sum(outer(w, w) * D) / 2
}

## ---------- 4. 系统发育多样性（probability-weighted） ---------------------

# 概率加权 Faith's PD（Reese-Tucker 2024 修订版本简版）：
#   PD = sum_{e in tree.edges} edge_length_e * (1 - prod_{sp under e}(1 - psi_sp))
div_pd_prob <- function(psi, phy, species_match) {
  if (!requireNamespace("ape", quietly = TRUE)) return(NA_real_)
  if (is.null(phy) || is.null(species_match)) return(NA_real_)
  psi <- as.numeric(psi)
  if (all(is.na(psi))) return(NA_real_)
  # species_match: 长度=length(psi)，元素是 phy$tip.label 的索引（NA 表示无匹配）
  m <- species_match
  edge_lengths <- phy$edge.length
  edge_descend <- phangorn::Descendants(phy, phy$edge[, 2], type = "tips")
  contrib <- vapply(seq_along(edge_lengths), function(i) {
    tips <- edge_descend[[i]]
    sp_idx <- which(!is.na(m) & m %in% tips)
    if (length(sp_idx) == 0) return(0)
    1 - prod(1 - pmin(pmax(psi[sp_idx], 0), 1))
  }, numeric(1))
  sum(edge_lengths * contrib)
}

# probability-weighted MPD（基于 cophenetic 距离）
div_mpd_prob <- function(psi, phy_dist_mat, species_match, eps = 1e-12) {
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps & !is.na(species_match)
  if (sum(ok) < 2) return(NA_real_)
  m   <- species_match[ok]
  w   <- psi[ok]; w <- w / sum(w)
  D   <- phy_dist_mat[m, m, drop = FALSE]
  diag(D) <- NA
  sum(outer(w, w) * D, na.rm = TRUE) / sum(outer(w, w) * !is.na(D))
}

## ---------- 5. 一站式：一个 grid × period 的 psi 向量 → 一行指标 ----------

community_metrics_one <- function(psi, trait_mat, dist_mat,
                                   phy = NULL, phy_dist = NULL,
                                   species_match = NULL) {
  tax <- div_taxonomic(psi)
  fn  <- div_functional(psi, trait_mat)
  raq <- div_rao_q(psi, dist_mat)
  pdv <- div_pd_prob(psi, phy, species_match)
  mpv <- div_mpd_prob(psi, phy_dist, species_match)
  c(tax,
    fn$cwm, fn$cwsd,
    trait_volume = unname(fn$trait_volume),
    trait_dispersion = unname(fn$trait_dispersion),
    rao_q = unname(raq),
    pd_prob = unname(pdv),
    mpd_prob = unname(mpv))
}

## ---------- 6. Temporal beta：Baselga 分解（基于 abundance/probability） --
##  连续 period 之间 (psi_t-1, psi_t)。
##  返回：bray, soerensen, turnover (Baselga), nestedness (Baselga)
div_temporal_beta_pair <- function(p_prev, p_curr) {
  p_prev <- as.numeric(p_prev); p_curr <- as.numeric(p_curr)
  ok <- is.finite(p_prev) & is.finite(p_curr)
  p_prev <- p_prev[ok]; p_curr <- p_curr[ok]
  if (length(p_prev) == 0) return(c(bray = NA, soerensen = NA,
                                     turnover = NA, nestedness = NA,
                                     delta_richness = NA))
  A <- pmin(p_prev, p_curr)              # shared
  B <- pmax(p_prev - p_curr, 0)          # lost
  C <- pmax(p_curr - p_prev, 0)          # gained
  bray <- (sum(B) + sum(C)) / (sum(p_prev) + sum(p_curr))
  soer <- 1 - 2 * sum(A) / (sum(p_prev) + sum(p_curr))
  bc <- min(sum(B), sum(C))
  turn <- if ((sum(A) + bc) > 0)
    bc / (sum(A) + bc) else NA_real_
  nest <- soer - turn
  c(bray = bray, soerensen = soer,
    turnover = turn, nestedness = nest,
    delta_richness = sum(p_curr) - sum(p_prev))
}

## ============================================================================
## v3 扩展：fundiversity FD + picante PSV/NRI/NTI + 自写 FSpe / FOri
## ============================================================================

## ---- FD via fundiversity（基于 trait PCoA 坐标 + abundance/probability 权重）

div_fd_extended <- function(psi, trait_pcoa_coords, trait_dist, eps = 1e-6) {
  if (!requireNamespace("fundiversity", quietly = TRUE)) {
    return(c(fric = NA_real_, feve = NA_real_, fdiv = NA_real_,
             fdis = NA_real_, fraoq = NA_real_,
             fspe = NA_real_, fori = NA_real_))
  }
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps
  if (sum(ok) < 4) return(c(fric=NA, feve=NA, fdiv=NA, fdis=NA, fraoq=NA, fspe=NA, fori=NA))
  sp_in <- rownames(trait_pcoa_coords)[ok]
  coords <- trait_pcoa_coords[ok, , drop = FALSE]
  w_mat  <- matrix(psi[ok], nrow = 1, dimnames = list("comm", sp_in))

  fric <- tryCatch(fundiversity::fd_fric(coords, w_mat)$FRic, error = function(e) NA_real_)
  feve <- tryCatch(fundiversity::fd_feve(coords, w_mat)$FEve, error = function(e) NA_real_)
  fdiv <- tryCatch(fundiversity::fd_fdiv(coords, w_mat)$FDiv, error = function(e) NA_real_)
  fdis <- tryCatch(fundiversity::fd_fdis(coords, w_mat)$FDis, error = function(e) NA_real_)
  fraoq <- tryCatch({
    D <- as.matrix(trait_dist)[ok, ok, drop = FALSE]
    w <- psi[ok]; w <- w / sum(w)
    sum(outer(w, w) * D) / 2
  }, error = function(e) NA_real_)

  # FSpe: 群落到性状空间整体重心的距离（特化度）
  centroid_all <- colMeans(trait_pcoa_coords)
  w <- psi[ok]; w <- w / sum(w)
  fspe <- sum(w * sqrt(rowSums((t(t(coords) - centroid_all))^2)))

  # FOri: 群落内每物种到其最近邻的距离的加权均值（原创度）
  if (sum(ok) >= 2) {
    D_in <- as.matrix(stats::dist(coords))
    diag(D_in) <- Inf
    nearest <- apply(D_in, 1, min)
    fori <- sum(w * nearest)
  } else fori <- NA_real_

  c(fric = as.numeric(fric), feve = as.numeric(feve), fdiv = as.numeric(fdiv),
    fdis = as.numeric(fdis), fraoq = as.numeric(fraoq),
    fspe = as.numeric(fspe), fori = as.numeric(fori))
}

## ---- 系统发育扩展：PSV / PSR / NRI / NTI（基于 picante）

div_phylo_extended <- function(psi, phy, species_match, phy_dist, eps = 1e-6) {
  if (is.null(phy) || is.null(species_match) || is.null(phy_dist)) {
    return(c(psv = NA, psr = NA, mpd_ses = NA, mntd_ses = NA,
             nri = NA, nti = NA))
  }
  psi <- as.numeric(psi)
  ok  <- is.finite(psi) & psi > eps & !is.na(species_match)
  if (sum(ok) < 2) return(c(psv=NA, psr=NA, mpd_ses=NA, mntd_ses=NA, nri=NA, nti=NA))

  # PSV / SES.MPD / SES.MNTD：picante 要求 community 矩阵 ≥ 2 行。
  # 用 2 行（focal community + 全员 1 的 reference）作为 workaround，取 row 1。
  m <- species_match[ok]; tips <- phy$tip.label[m]
  pres <- matrix(0, 2, ape::Ntip(phy),
                 dimnames = list(c("comm","reference"), phy$tip.label))
  pres[1, tips] <- psi[ok]
  pres[2, ]     <- 1                # 参考全员存在

  psv_v <- tryCatch(picante::psv(pres, phy, compute.var = FALSE)$PSVs[1],
                    error = function(e) NA_real_)
  psr_v <- tryCatch(picante::psr(pres, phy, compute.var = FALSE)$PSR[1],
                    error = function(e) NA_real_)

  ses_mpd <- tryCatch({
    s <- picante::ses.mpd(pres, phy_dist, null.model = "taxa.labels",
                           abundance.weighted = TRUE, runs = 19)
    list(mpd_ses = s$mpd.obs.z[1], nri = -s$mpd.obs.z[1])
  }, error = function(e) list(mpd_ses = NA, nri = NA))

  ses_mntd <- tryCatch({
    s <- picante::ses.mntd(pres, phy_dist, null.model = "taxa.labels",
                            abundance.weighted = TRUE, runs = 19)
    list(mntd_ses = s$mntd.obs.z[1], nti = -s$mntd.obs.z[1])
  }, error = function(e) list(mntd_ses = NA, nti = NA))

  c(psv = as.numeric(psv_v), psr = as.numeric(psr_v),
    mpd_ses = as.numeric(ses_mpd$mpd_ses),
    mntd_ses = as.numeric(ses_mntd$mntd_ses),
    nri = as.numeric(ses_mpd$nri),
    nti = as.numeric(ses_mntd$nti))
}

## ---- Phylo-β（PhyloSor）相邻 period 间

div_phylosor_pair <- function(p_prev, p_curr, phy, species_match, eps = 1e-6) {
  if (is.null(phy) || is.null(species_match)) return(NA_real_)
  m1 <- which(p_prev > eps & !is.na(species_match))
  m2 <- which(p_curr > eps & !is.na(species_match))
  if (length(m1) < 2 || length(m2) < 2) return(NA_real_)
  pres <- matrix(0, 2, ape::Ntip(phy),
                 dimnames = list(c("t1","t2"), phy$tip.label))
  pres[1, phy$tip.label[species_match[m1]]] <- 1
  pres[2, phy$tip.label[species_match[m2]]] <- 1
  v <- tryCatch(picante::phylosor(pres, phy), error = function(e) NA_real_)
  if (length(v) == 0) NA_real_ else as.numeric(v[1])
}
