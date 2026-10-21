#!/usr/bin/env Rscript
# ============================================================
# 31_protected_area_effectiveness.R
#
# Scientific question / 科学问题:
#   中国自然保护区是否延缓了鸟类群落的功能同质化？
#   保护区选址高度非随机（多在高海拔、偏远、低压力区，Joppa & Pfaff 2009），
#   因此"区内 vs 区外"的朴素比较必然有偏。本脚本用倾向性得分匹配 +
#   双重差分，估计保护区对多维生物多样性趋势的因果效应，并做保护空缺分析。
#   Do China's nature reserves slow avian functional homogenization?
#
# Objective / 分析目标:
#   1) 网格保护区覆盖率（按级别/类型/建立年份分层）
#   2) 匹配设计：以基线协变量匹配区内外网格，比较多维趋势（ATT）
#   3) 双重差分：2000-2024 年新建保护区（n≈285）的建立前后对比
#   4) 类型分解：湿地类保护区 vs 森林类保护区的差异效应（检验水鸟假说）
#   5) 保护空缺：物种代表性、功能独特性空缺、同质化热点空缺
#
# Input / 输入:
#   external: protected_areas/全国自然保护区名录+矢量边界/保护区.shp (n=1028, EPSG:4326)
#   derived : grid_environment<GRID_TAG>_v3.rds  （匹配协变量：海拔、基线HFI、气候、地被）
#   results : table_trend_summary_*  / table_diversity_summary_*  / table_cti_*（脚本30产出）
#   derived : psi_samples_thinned_*  （物种代表性与功能空缺）
# Output / 输出:
#   results: table_pa_grid_coverage_<run_label>.csv
#   results: table_pa_matching_att_<run_label>.csv         匹配后处理效应
#   results: table_pa_did_<run_label>.csv                  双重差分
#   results: table_pa_type_effects_<run_label>.csv         按保护类型
#   results: table_conservation_gap_<run_label>.csv        物种/功能代表性空缺
#
# Key assumptions / 关键假设:
#   - 可忽略性：匹配变量已捕捉保护区选址的主要非随机来源（海拔、基线人类压力、气候、地被）
#     —— 无法排除未观测混杂，故结论表述为"与保护相关的差异"而非严格因果
#   - 覆盖率按面积加权；保护区矢量仅 1028/3376 条（以国家级与大面积为主），
#     故结果偏向大型/高级别保护区，须在文中声明
# Packages: sf, dplyr, MatchIt, readr
# References: Andam et al. 2008 PNAS; Joppa & Pfaff 2009 PLoS ONE;
#             Ferraro & Hanauer 2014 Annu Rev Environ Resour; Geldmann et al. 2019 PNAS
# ============================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))

ensure_v3_dirs()
is_pilot  <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
log_time("31", "Protected-area effectiveness and conservation gaps")

# 优先读 UTF-8 的 GeoPackage（跨平台安全）；缺失时回退到原始 shapefile。
# 原始 shp 为 GBK 编码，若不显式声明，Linux(UTF-8 locale) 下中文字段会乱码，
# 导致属性连接静默失败。Prefer UTF-8 GeoPackage; the source shapefile is GBK-encoded.
PA_DIR <- file.path(DIRS$external, "protected_areas")
PA_GPKG <- file.path(PA_DIR, "china_nature_reserves_utf8.gpkg")
PA_SHP  <- file.path(PA_DIR, "全国自然保护区名录+矢量边界", "保护区.shp")
if (!file.exists(PA_GPKG) && !file.exists(PA_SHP))
  stop("[31] PA boundaries not found. Expected:\n  ", PA_GPKG, "\n  or ", PA_SHP)

# ── 1. 网格 × 保护区覆盖率 ───────────────────────────────────────────
pa <- if (file.exists(PA_GPKG)) {
  message("[31] reading PA boundaries from GeoPackage (UTF-8)")
  st_read(PA_GPKG, quiet = TRUE)
} else {
  message("[31] reading PA boundaries from shapefile with explicit GBK encoding")
  st_read(PA_SHP, quiet = TRUE, options = "ENCODING=GBK")
}
pa <- st_make_valid(pa)
# 编码自检：保护区名称若出现替换字符则说明编码仍不对，宁可停下也不要静默出错
nm_col <- intersect(c("保护区名称", "name"), names(pa))[1]
if (!is.na(nm_col) && any(grepl("�", pa[[nm_col]]), na.rm = TRUE))
  stop("[31] PA name field is mojibake — encoding is wrong; use the UTF-8 GeoPackage.")
pa$year <- suppressWarnings(as.numeric(substr(as.character(pa$`始建时间`), 1, 4)))
message(sprintf("[31] PAs: %d polygons; pre-2000 %d; 2000-2024 %d",
                nrow(pa), sum(pa$year < 2000, na.rm = TRUE),
                sum(pa$year >= 2000 & pa$year <= 2024, na.rm = TRUE)))

# 载入与建模完全一致的网格 sf（02/03 脚本生成并落盘）
# Load the same grid sf used for modelling (written by scripts 02/03).
grid_sf <- safe_read(v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds"))
if (is.null(grid_sf)) grid_sf <- safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"))
if (is.null(grid_sf)) stop("[31] china_grid rds not found; run 02 first.")
grid_sf <- st_transform(grid_sf, st_crs(pa))
# 仅保留进入模型的网格，保证与结果表逐格对齐 / keep modelled grids only
surv <- safe_read(v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds"))
if (!is.null(surv) && !is.null(surv$sites))
  grid_sf <- grid_sf[grid_sf$grid_cell %in% surv$sites, ]
message(sprintf("[31] grid cells used: %d", nrow(grid_sf)))

# 逐网格计算被保护区覆盖的面积比例（等积投影下求面积，避免经纬度面积失真）
# Area fraction of each grid covered by PAs, computed in an equal-area projection.
cover_frac <- function(pa_subset) {
  n <- nrow(grid_sf)
  if (nrow(pa_subset) == 0) return(rep(0, n))
  u <- st_union(st_make_valid(st_geometry(pa_subset)))
  g_ea <- project_china_albers(grid_sf)          # utils_spatial 提供的等积投影
  u_ea <- st_transform(st_sf(geometry = u), st_crs(g_ea))
  grid_area <- as.numeric(st_area(g_ea))
  hit <- st_intersects(g_ea, u_ea, sparse = FALSE)[, 1]
  out <- rep(0, n)
  if (any(hit)) {
    inter <- st_intersection(st_geometry(g_ea)[hit], st_geometry(u_ea))
    out[hit] <- as.numeric(st_area(inter))
  }
  pmin(out / grid_area, 1)
}

cov_tbl <- tibble(
  grid_cell   = grid_sf$grid_cell,
  pa_frac_all = cover_frac(pa),
  pa_frac_pre2000  = cover_frac(pa[which(pa$year <  2000), ]),
  pa_frac_post2000 = cover_frac(pa[which(pa$year >= 2000), ]),
  pa_frac_national = cover_frac(pa[which(pa$`级别` == "国家级"), ]),
  pa_frac_wetland  = cover_frac(pa[which(pa$`类型` %in% c("内陆湿地", "海洋海岸")), ]),
  pa_frac_forest   = cover_frac(pa[which(pa$`类型` == "森林生态"), ]))
# 规范列名：下游脚本（33/34）以 pa_frac 为准，pa_frac_all 保留兼容旧结果
cov_tbl$pa_frac <- cov_tbl$pa_frac_all

# ── 1b. 网格的保护区设立年份（事件研究 / 动态 DiD 的处理时点）──────────
# 逐网格取"与之相交的保护区中最早的设立年份"。
#   pa_year_min       : 任何相交保护区的最早年份（哪怕只压到一角）
#   pa_year_first_major: 仅计入贡献网格面积 >= MAJOR_FRAC 的保护区，
#                        避免一小块老保护区的边角把处理时点错误地提前——
#                        事件研究应以此列为准。
# Earliest reserve establishment year per grid; the "major" variant ignores
# slivers so the event-study treatment date is not set by a marginal overlap.
MAJOR_FRAC <- 0.01
message("[31] 计算网格保护区设立年份 ...")
g_ea   <- project_china_albers(grid_sf)
pa_ea  <- st_transform(pa, st_crs(g_ea))
grid_a <- as.numeric(st_area(g_ea))
hits   <- st_intersects(g_ea, pa_ea)          # 稀疏列表：每网格相交的保护区索引

yr_min <- rep(NA_real_, nrow(grid_sf))
yr_maj <- rep(NA_real_, nrow(grid_sf))
n_pa   <- lengths(hits)
for (i in seq_len(nrow(grid_sf))) {
  idx <- hits[[i]]
  if (!length(idx)) next
  yrs <- pa_ea$year[idx]
  if (all(is.na(yrs))) next
  yr_min[i] <- suppressWarnings(min(yrs, na.rm = TRUE))
  # 各相交保护区对该网格的面积贡献
  inter <- suppressWarnings(
    st_intersection(st_geometry(g_ea)[i], st_geometry(pa_ea)[idx]))
  if (length(inter) == 0) next
  fr <- as.numeric(st_area(inter)) / grid_a[i]
  keep <- which(fr >= MAJOR_FRAC)
  if (length(keep)) {
    yk <- yrs[seq_along(fr)][keep]
    if (any(!is.na(yk))) yr_maj[i] <- suppressWarnings(min(yk, na.rm = TRUE))
  }
}
yr_min[!is.finite(yr_min)] <- NA_real_
yr_maj[!is.finite(yr_maj)] <- NA_real_
cov_tbl$n_pa_overlapping  <- n_pa
cov_tbl$pa_year_min       <- yr_min
cov_tbl$pa_year_first_major <- yr_maj
# 事件研究可用的处理组：2000-2024 期间首次获得实质保护的网格
cov_tbl$newly_protected_2000_2024 <-
  !is.na(yr_maj) & yr_maj >= 2000 & yr_maj <= 2024

write_csv(cov_tbl, paste0(v3_file("results", paste0("table_pa_grid_coverage_", run_label)), ".csv"))
message(sprintf("[31] 网格平均保护区覆盖率 %.1f%%；覆盖率>10%% 的网格 %d 个",
                100 * mean(cov_tbl$pa_frac_all), sum(cov_tbl$pa_frac_all > 0.10)))
message(sprintf("[31] 有设立年份的网格 %d 个；2000-2024 期间新获实质保护 %d 个（事件研究处理组）",
                sum(!is.na(cov_tbl$pa_year_first_major)),
                sum(cov_tbl$newly_protected_2000_2024)))
if (sum(cov_tbl$newly_protected_2000_2024) < 20)
  message("[31] ⚠ 事件研究处理组样本偏少，动态 DiD 结果需谨慎解读")

# ── 2. 组装分析数据：覆盖率 + 匹配协变量 + 多维趋势结果 ──────────────
genv <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(genv)) genv <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))

trend <- read_csv(paste0(v3_file("results",
          paste0("table_trend_summary_", run_label)), "_extended.csv"), show_col_types = FALSE)
wide <- trend |>
  filter(method %in% c("theil_sen", "ols")) |>
  select(grid_cell, metric, mean) |>
  tidyr::pivot_wider(names_from = metric, values_from = mean, names_prefix = "trend_")

dat <- cov_tbl |>
  left_join(genv, by = "grid_cell") |>
  left_join(wide, by = "grid_cell") |>
  mutate(treated = as.integer(pa_frac_all >= 0.10))   # 处理组定义：覆盖率≥10%

# ── 3. 倾向性得分匹配：控制保护区选址偏倚 ────────────────────────────
# 协变量选择依据 Joppa & Pfaff (2009)：保护区偏好高海拔、陡峭、偏远、低压力区
COVS <- c("elev_mean", "elev_sd", "bio1", "bio12", "hfi_2000",
          "landcover_trees", "landcover_cropland", "landcover_built",
          "npp_mean", "centroid_lon", "centroid_lat")
COVS <- intersect(COVS, names(dat))
mdat <- dat |> filter(if_all(all_of(c(COVS, "treated")), ~ !is.na(.)))
message(sprintf("[31] 匹配样本：处理组 %d，对照组 %d",
                sum(mdat$treated == 1), sum(mdat$treated == 0)))

if (!requireNamespace("MatchIt", quietly = TRUE)) {
  message("[31] MatchIt 未安装，跳过匹配；install.packages('MatchIt')")
} else {
  f <- as.formula(paste("treated ~", paste(COVS, collapse = " + ")))
  m <- MatchIt::matchit(f, data = mdat, method = "nearest",
                        distance = "glm", caliper = 0.2, ratio = 1)
  print(summary(m)$sum.matched[, c("Means Treated", "Means Control", "Std. Mean Diff.")])

  # 导出匹配前后的协变量平衡（|SMD|<0.1 视为平衡）。这是准实验设计的
  # 必要证据：不给平衡表，审稿人无法判断 ATT 是否可信。
  # Export covariate balance before/after matching; |SMD| < 0.1 = balanced.
  sm <- summary(m)
  bal <- data.frame(
    variable  = rownames(sm$sum.all),
    smd_before = sm$sum.all[, "Std. Mean Diff."],
    smd_after  = sm$sum.matched[match(rownames(sm$sum.all),
                                      rownames(sm$sum.matched)), "Std. Mean Diff."],
    stringsAsFactors = FALSE)
  bal$balanced_after <- abs(bal$smd_after) < 0.1
  write_csv(bal, paste0(v3_file("results", paste0("table_pa_matching_balance_", run_label)), ".csv"))
  message(sprintf("[31] 匹配平衡：%d/%d 个协变量在匹配后 |SMD|<0.1",
                  sum(bal$balanced_after, na.rm = TRUE), nrow(bal)))
  if (any(!bal$balanced_after, na.rm = TRUE))
    message("[31] ⚠ 未平衡的协变量：",
            paste(bal$variable[!bal$balanced_after & !is.na(bal$balanced_after)], collapse = ", "),
            " —— ATT 解释需谨慎")

  md <- MatchIt::match.data(m)

  OUTCOMES <- grep("^trend_", names(md), value = TRUE)
  att <- do.call(rbind, lapply(OUTCOMES, function(o) {
    if (all(is.na(md[[o]]))) return(NULL)
    tt <- t.test(md[[o]][md$treated == 1], md[[o]][md$treated == 0])
    data.frame(outcome = o,
               mean_protected   = unname(tt$estimate[1]),
               mean_unprotected = unname(tt$estimate[2]),
               ATT = unname(tt$estimate[1] - tt$estimate[2]),
               ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
               p_value = tt$p.value, n_matched = nrow(md),
               stringsAsFactors = FALSE)
  }))
  print(att)
  write_csv(att, paste0(v3_file("results", paste0("table_pa_matching_att_", run_label)), ".csv"))
  message("[31] >>> 若 trend_fric_prob / trend_rao_q 的 ATT 显著为正，")
  message("[31] >>> 说明保护区延缓了功能同质化——这是保护成效的核心证据。")

  # ── 4. 按保护类型分解：湿地 vs 森林（检验水鸟驱动假说）──────────────
  ty <- do.call(rbind, lapply(c("pa_frac_wetland", "pa_frac_forest"), function(v) {
    d2 <- dat |> mutate(tr = as.integer(.data[[v]] >= 0.10)) |>
      filter(if_all(all_of(COVS), ~ !is.na(.)))
    if (sum(d2$tr) < 20) return(NULL)
    f2 <- as.formula(paste("tr ~", paste(COVS, collapse = " + ")))
    m2 <- MatchIt::matchit(f2, data = d2, method = "nearest", distance = "glm",
                           caliper = 0.2, ratio = 1)
    d3 <- MatchIt::match.data(m2)
    do.call(rbind, lapply(grep("^trend_", names(d3), value = TRUE), function(o) {
      if (all(is.na(d3[[o]]))) return(NULL)
      tt <- t.test(d3[[o]][d3$tr == 1], d3[[o]][d3$tr == 0])
      data.frame(pa_type = v, outcome = o, ATT = unname(diff(rev(tt$estimate))),
                 p_value = tt$p.value, n_treated = sum(d3$tr), stringsAsFactors = FALSE)
    }))
  }))
  print(ty); write_csv(ty, paste0(v3_file("results", paste0("table_pa_type_effects_", run_label)), ".csv"))
}

# ── 5. 双重差分：2000-2024 年新建保护区 ──────────────────────────────
# 处理组 = 研究期内新建且期前无保护；结果 = 各期校正物种数/功能体积
div <- read_csv(paste0(v3_file("results",
        paste0("table_diversity_summary_", run_label)), "_extended.csv"), show_col_types = FALSE)
did_dat <- div |>
  filter(metric %in% c("corrected_richness", "fric_prob", "rao_q")) |>
  select(grid_cell, period, metric, value = mean) |>
  left_join(cov_tbl, by = "grid_cell") |>
  mutate(treated_new = as.integer(pa_frac_post2000 >= 0.05 & pa_frac_pre2000 < 0.01),
         post = as.integer(period %in% c("P4", "P5")))
did <- do.call(rbind, lapply(unique(did_dat$metric), function(mm) {
  d <- filter(did_dat, metric == mm)
  if (sum(d$treated_new) < 20) return(NULL)
  fit <- lm(value ~ treated_new * post + pa_frac_pre2000, data = d)
  s <- summary(fit)$coefficients
  rn <- "treated_new:post"
  if (!rn %in% rownames(s)) return(NULL)
  data.frame(metric = mm, did_estimate = s[rn, 1], se = s[rn, 2],
             p_value = s[rn, 4], n_treated_grids = sum(d$treated_new) / 5,
             stringsAsFactors = FALSE)
}))
print(did); write_csv(did, paste0(v3_file("results", paste0("table_pa_did_", run_label)), ".csv"))
message("[31] >>> DiD 为新建保护区提供准实验证据；正系数=保护后多样性相对改善。")

# ── 6. 保护空缺分析：物种代表性 + 功能独特性 + 同质化热点 ────────────
psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
if (!is.null(psi_obj) && length(dim(psi_obj$psi_samples_thinned)) >= 4) {
  psi <- apply(psi_obj$psi_samples_thinned, c(2, 3, 4), mean)
  sp <- psi_obj$species; sites <- psi_obj$sites; np <- dim(psi)[3]
  pf <- cov_tbl$pa_frac_all[match(sites, cov_tbl$grid_cell)]; pf[is.na(pf)] <- 0

  # 物种代表性：占域加权的保护区覆盖比例（对照昆明-蒙特利尔 30% 目标）
  rep_tbl <- data.frame(
    species = sp,
    protected_fraction = sapply(seq_along(sp), function(s) {
      w <- psi[s, , np]; sum(w * pf) / sum(w) }),
    stringsAsFactors = FALSE)
  rep_tbl$meets_30pct <- rep_tbl$protected_fraction >= 0.30
  message(sprintf("[31] 物种代表性：中位保护比例 %.1f%%；达到 30%% 目标的物种 %d/%d",
                  100 * median(rep_tbl$protected_fraction),
                  sum(rep_tbl$meets_30pct), nrow(rep_tbl)))

  # 功能独特性空缺：功能越独特的物种是否保护越不足？
  tr_path <- paste0(v3_file("results", paste0("table_species_trend_traits_", run_label)), "_extended.csv")
  if (file.exists(tr_path)) {
    tr <- read_csv(tr_path, show_col_types = FALSE)
    if ("habitat_breadth" %in% names(tr)) {
      rep_tbl <- left_join(rep_tbl, tr |> select(species, habitat_breadth, trend_slope), by = "species")
      ct <- cor.test(rep_tbl$habitat_breadth, rep_tbl$protected_fraction,
                     method = "spearman", exact = FALSE)
      message(sprintf("[31] 栖息地宽度 vs 保护比例：rho = %.3f, p = %.3g", ct$estimate, ct$p.value))
      message("[31] >>> rho 显著为正 = 泛化种反而受更好保护 = 特化种保护空缺。")
    }
  }
  write_csv(rep_tbl, paste0(v3_file("results", paste0("table_conservation_gap_", run_label)), ".csv"))
} else message("[31] 无 4D psi，跳过代表性空缺分析。")

log_time("31", "DONE")
