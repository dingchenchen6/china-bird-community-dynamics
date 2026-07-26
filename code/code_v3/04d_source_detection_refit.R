#!/usr/bin/env Rscript
# ============================================================
# 04d_source_detection_refit.R
#
# Scientific question / 科学问题:
#   把数据来源(平台)异质性纳入检测子模型后,群落占域轨迹是否改变?
#   即:eBird/GBIF 与中国观鸟记录中心(CBRC)协议不同带来的检测差异,
#   是否被误当成占域变化? 用 ΔWAIC 与斜率变化回答。
#   Does adding a data-source detection term change the community trajectory?
#
# Objective / 分析目标:
#   在 04 的 stMsPGOcc 基础上,仅给 detection 子模型增加一个来源协变量
#   (每个 site×primary×secondary 单元的非-birdwatch 记录占比),其余设定完全一致,
#   重跑并用 spOccupancy::waicOcc 比较 base vs +source。
#
# Input / 输入:
#   derived: stMsPGOcc_fit_<run_label>.rds        (取 y / occ.covs / det.covs / coords / X)
#   derived: survey_history<GRID_TAG>_v3.rds       (取 replicate 级 source;见 >>> ADJUST)
# Output / 输出:
#   derived: stMsPGOcc_fit_<run_label>_srcdet.rds
#   results: table_source_detection_waic_<run_label>.csv  (base/+source 的 WAIC、ΔWAIC)
#
# Key assumption / 关键假设:
#   source 协变量在 detection 层可识别;occ.formula/priors/inits/tuning/MCMC
#   与 04 完全一致(否则 WAIC 不可比)。
# Packages: spOccupancy
#
# ⚠️ 三项 server 补算中最重的一项:需按你实际 survey_history 的字段名核对
#    两处 `## >>> ADJUST`。其余逻辑与 04 对齐。
# ============================================================

suppressPackageStartupMessages({
  library(spOccupancy)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project",
            "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

ensure_v3_dirs()
is_pilot  <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
set.seed(if (exists("SEED")) SEED else 20260602L)
log_time("04d", "Source-detection refit + WAIC comparison")

# ── 1. 复用 04 的 base fit 里保存的模型输入 ─────────────────────────
fit_base <- safe_read(v3_file("derived", paste0("stMsPGOcc_fit_", run_label), "rds"))
if (is.null(fit_base)) stop("[04d] base fit not found; run 04 first.")

y        <- fit_base$y            # [species, site, primary, secondary]
occ_covs <- fit_base$occ.covs     # list / data.frame of [site(,primary)]
det_covs <- fit_base$det.covs     # list of [site, primary, secondary]
coords   <- fit_base$coords
dy <- dim(y); n_sites <- dy[2]; n_primary <- dy[3]; n_secondary <- dy[4]
message(sprintf("[04d] y dims: species=%d sites=%d primary=%d secondary=%d",
                dy[1], n_sites, n_primary, n_secondary))

# occ.formula:与 04 一致(reformulate(c(occ_vars, "year_scaled")))。
# 从 base fit 的设计矩阵列名重建,去掉截距,避免手误。
occ_terms <- setdiff(colnames(fit_base$X), c("(Intercept)"))
occ_formula <- reformulate(occ_terms)
det_formula_src <- ~ log_events + log_duration + has_duration + src_prop

# ── 2. 构造来源协变量 src_prop[site, primary, secondary] ──────────────
#   定义:该检测单元中"非 birdwatch(即 eBird/GBIF)记录占比"(0=全 CBRC, 1=全外源)。
survey <- safe_read(v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds"))
if (is.null(survey)) survey <- safe_read(v3_file("derived", "survey_history_v3", "rds"))
if (is.null(survey)) stop("[04d] survey_history not found.")

## >>> ADJUST #1: 取出 replicate 级长表,需含列:
##     site_index, block_id(=primary), year_in_block(=secondary), source
##     04 第139行即用 rr$site_index/rr$block_id/rr$year_in_block 填 det 数组,
##     若该长表在 survey 里叫别的名字,请改这一行。
recs <- survey$records
if (is.null(recs)) recs <- survey$long
stopifnot(all(c("site_index", "block_id", "year_in_block") %in% names(recs)))
if (!"source" %in% names(recs))
  stop("[04d] 'source' 列不在 survey records 中。需在 02_build_survey_history.R \n",
       "        为每条 detection 记录保留 source(birdwatch/ebird/gbif),或提供含 source 的原始 events 映射。")

# 聚合:每个 (site,primary,secondary) 的非-birdwatch 占比
recs$is_ext <- as.integer(tolower(recs$source) != "birdwatch")   # ebird/gbif = 1
agg <- aggregate(is_ext ~ site_index + block_id + year_in_block, data = recs, FUN = mean)

src_arr <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
for (i in seq_len(nrow(agg))) {
  si <- agg$site_index[i]; bi <- agg$block_id[i]; yi <- agg$year_in_block[i]
  if (si >= 1 && si <= n_sites && bi >= 1 && bi <= n_primary &&
      yi >= 1 && yi <= n_secondary)
    src_arr[si, bi, yi] <- agg$is_ext[i]
}
# 缺失单元:用已访问单元的均值填补,并标准化(与 04 检测协变量同风格)
visited <- is.finite(det_covs$log_events)
src_arr[is.na(src_arr) & visited] <- mean(src_arr, na.rm = TRUE)
src_arr[!is.finite(src_arr)] <- mean(src_arr, na.rm = TRUE)
src_arr <- (src_arr - mean(src_arr, na.rm = TRUE)) / sd(as.vector(src_arr), na.rm = TRUE)

det_covs$src_prop <- src_arr
message(sprintf("[04d] src_prop built (scaled); non-NA cells: %d", sum(is.finite(src_arr))))

# ── 3. 组装 data,重跑(设定必须与 04 一致)─────────────────────────────
data_list <- list(y = y, occ.covs = occ_covs, det.covs = det_covs, coords = coords)
z_init <- apply(y, c(1, 2, 3), function(a) as.integer(any(a == 1, na.rm = TRUE)))

## >>> ADJUST #2: 以下 inits / priors / tuning / n.batch 必须与 04_run_stMsPGOcc_main.R
##     第 357–405 行完全一致(否则 WAIC 不可比)。建议直接从 04 复制粘贴。
n_batch      <- if (is_pilot) 200L else 1600L   # 04: n_batch × 25 = 总迭代;请对齐 04
n_burn       <- FULL_N_BURN
n_thin       <- FULL_N_THIN
n_chains     <- FULL_N_CHAINS
n_omp        <- as.integer(Sys.getenv("V3_OMP_THREADS", "1"))
COV_MODEL    <- if (exists("COV_MODEL")) COV_MODEL else "exponential"
N_NEIGHBORS  <- if (exists("N_NEIGHBORS")) N_NEIGHBORS else 15L

inits  <- list(z = z_init)                       # >>> 复制 04 的完整 inits
priors <- list()                                  # >>> 复制 04 的完整 priors
tuning <- list(phi = 1)                           # >>> 复制 04 的完整 tuning
if (length(priors) == 0L)
  warning("[04d] priors 为空占位——请从 04 复制真实 priors 后再跑,否则结果不可比。")

fit_src <- stMsPGOcc(
  occ.formula   = occ_formula,
  det.formula   = det_formula_src,
  data          = data_list,
  inits         = inits,
  priors        = priors,
  tuning        = tuning,
  cov.model     = COV_MODEL,
  NNGP          = TRUE,
  n.neighbors   = N_NEIGHBORS,
  n.batch       = n_batch,
  batch.length  = 25,
  n.burn        = n_burn,
  n.thin        = n_thin,
  n.chains      = n_chains,
  ar1           = TRUE,
  n.omp.threads = n_omp,
  verbose       = TRUE)

saveRDS(fit_src, v3_file("derived", paste0("stMsPGOcc_fit_", run_label, "_srcdet"), "rds"),
        compress = "xz")

# ── 4. WAIC 比较 ─────────────────────────────────────────────────────
w_base <- as.numeric(waicOcc(fit_base)["WAIC"])
w_src  <- as.numeric(waicOcc(fit_src)["WAIC"])
res <- data.frame(
  model = c("base (no source)", "+ source detection"),
  WAIC  = c(w_base, w_src),
  delta_WAIC_vs_base = c(0, w_src - w_base))
print(res, row.names = FALSE)
write.csv(res, paste0(v3_file("results",
          paste0("table_source_detection_waic_", run_label)), ".csv"), row.names = FALSE)

message(sprintf("[04d] ΔWAIC (+source − base) = %.1f  (负值=source 项改善拟合)", w_src - w_base))
log_time("04d", "DONE")
