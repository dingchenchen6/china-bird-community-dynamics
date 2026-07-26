#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   v4 全局配置：在 v3 基础上修正版本沿袭、路径本地化、bug 防御。
#   所有 v4 脚本 source 此文件获取路径与参数，绝不执行计算。
#
# Objective / 分析目标:
#   单点定义 PROJECT_ROOT / DIRS / 模型参数 / 性状变换分类 / Nature 规范。
#
# Input data / 输入数据:
#   环境变量 BIRD_PROJECT_ROOT（可选，否则自动检测）
#
# Main workflow / 主要流程:
#   1. 定位项目根目录
#   2. 定义 v4 独立目录布局（不覆盖 v3）
#   3. 时间 / 空间 / MCMC / brms / 性状参数
#   4. 外部栅格 / 性状 / 系统发育路径（本地化优先）
#   5. Nature 出版图规范
#
# Key assumptions / 关键假设:
#   - v4 输出全部隔离到 code_v4/results_v4/figures_v4/data/derived_v4
#   - 如本地化文件缺失，自动 fallback 到 v3 derived 或邻近项目
#
# Main packages / 主要包:
#   无（仅 base R 路径操作）
#
# Output directory / 输出路径:
#   不产出文件；定义全局常量
# ============================================================

# ── 项目根目录（自动检测） ───────────────────────────────────────────
PROJECT_ROOT <- Sys.getenv(
  "BIRD_PROJECT_ROOT",
  if (dir.exists(file.path("~", "Documents", "New project",
                           "bird_dynamic_occupancy_analysis"))) {
    file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
  } else {
    getwd()
  }
)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)

# ── 目录布局（v4 全独立，v3/v2 仅作 fallback 数据源） ────────────────
DIRS <- list(
  code         = file.path(PROJECT_ROOT, "code_v4"),
  data_raw     = file.path(PROJECT_ROOT, "data"),
  derived      = file.path(PROJECT_ROOT, "data", "derived_v4"),
  results      = file.path(PROJECT_ROOT, "results_v4"),
  figures      = file.path(PROJECT_ROOT, "figures_v4"),
  logs         = file.path(PROJECT_ROOT, "logs_v4"),
  # fallback sources（read-only）
  v3_derived   = file.path(PROJECT_ROOT, "data", "derived_v3"),
  v3_results   = file.path(PROJECT_ROOT, "results_v3"),
  v2_derived   = file.path(PROJECT_ROOT, "data", "derived_v2"),
  v2_results   = file.path(PROJECT_ROOT, "results_v2"),
  external     = file.path(PROJECT_ROOT, "data", "external")
)

# ── 时间范围 ──────────────────────────────────────────────────────────
ANALYSIS_YR_LO   <- 2000L
ANALYSIS_YR_HI   <- 2024L
PERIOD_LENGTH    <- 5L          # 5 年一个 primary period
BREEDING_MONTHS  <- 4:8         # 4–8 月，可通过 env 覆盖

# ── 空间网格 ──────────────────────────────────────────────────────────
GRID_SIZE_KM     <- as.integer(Sys.getenv("V4_GRID_SIZE_KM", "100"))
if (length(GRID_SIZE_KM) == 0 || is.na(GRID_SIZE_KM)) GRID_SIZE_KM <- 100L
MIN_VISITS       <- as.integer(Sys.getenv("V4_MIN_VISITS", "3"))
MIN_SPECIES_GRID <- as.integer(Sys.getenv("V4_MIN_SPECIES_GRID", "5"))

# ── 运行标签 ──────────────────────────────────────────────────────────
GRID_TAG <- if (GRID_SIZE_KM != 100L) paste0("_", GRID_SIZE_KM, "km") else ""
RUN_LABEL <- Sys.getenv("V4_RUN_LABEL",
  paste0("v4_full_200sp_ar1_spatial", GRID_TAG))
PILOT_LABEL <- Sys.getenv("V4_PILOT_LABEL",
  paste0("v4_pilot_60sp_ar1_spatial", GRID_TAG))

# ── stMsPGOcc 模型参数 ─────────────────────────────────────────────
# Pilot（快速验证：本机 24 GB 可跑）
PILOT_N_BATCH   <- 225L
PILOT_N_BURN    <- 2000L
PILOT_N_THIN    <- 1L
PILOT_N_CHAINS  <- 4L
PILOT_MAX_SP    <- 20L

# Full（正式运行，建议在服务器 ≥ 256 GB）
FULL_N_BATCH   <- 400L
FULL_N_BURN    <- 5000L
FULL_N_THIN    <- 2L
FULL_N_CHAINS  <- 4L

# 空间参数（stMsPGOcc / NNGP）
N_NEIGHBORS     <- 5L
COV_MODEL       <- "exponential"   # "exponential", "spherical", "gaussian", "matern"

# MCMC 收敛阈值
RHAT_THRESHOLD  <- 1.05
ESS_THRESHOLD   <- 200L

# ── 后处理参数 ──────────────────────────────────────────────────────
PSI_MAX_DRAWS    <- 400L
POST_DRAWS_USE   <- 300L
TREND_DRAWS      <- 400L

# ── brms 参数 ──────────────────────────────────────────────────────
BRMS_ITER        <- 4000L
BRMS_WARMUP      <- 2000L
BRMS_CHAINS      <- 4L
BRMS_ADAPT_DELTA <- 0.99
BRMS_MAX_TREED   <- 15L
BRMS_SEED        <- 2024L

# v4 新增：horseshoe 稀疏先验参数（C6 修复）
BRMS_HS_DF       <- 1L
BRMS_HS_PAR      <- 0.1

# brms 空间 GP 基函数 k 三档敏感性
BRMS_GP_K_VEC    <- c(10L, 20L, 50L)

# ── 性状参数 ──────────────────────────────────────────────────────
TRAIT_VARS_BASIC <- c(
  "body_mass_g", "clutch_size", "longevity_y",
  "maturity_y", "avonet_hwi", "avonet_range_size"
)

TRAIT_VARS_EXTENDED <- c(
  "diet_specialization",
  "habitat_breadth"
)

TRAIT_VARS_REGRESSION_Z <- c(
  "z_body_mass", "z_hwi", "z_range_size", "z_clutch_size",
  "z_diet_specialization", "z_habitat_breadth"
)

TRAIT_VARS_ALL <- c(TRAIT_VARS_BASIC, TRAIT_VARS_EXTENDED)

# 性状变换分类（继承 v3 修复 #13）
TRAIT_VARS_LOG10  <- c("body_mass_g", "clutch_size", "longevity_y",
                       "maturity_y", "avonet_range_size")
TRAIT_VARS_NO_LOG <- c("avonet_hwi", "diet_specialization", "habitat_breadth")

# ── 多样性计算阈值（v4 新增：eps 敏感性） ─────────────────────────────
# eps 在 div_functional / div_phylogenetic_prob 中作为"存在概率"阈值
# v4 默认 0.05（相当于物种在该网格的占有率 ≥ 5% 才纳入加权）
# 17_sensitivity_eps_threshold.R 会扫描 1e-12 / 1e-6 / 0.01 / 0.05 / 0.10
PSI_EPS_DEFAULT  <- 0.05

# ── 外部数据路径（本地化优先，邻近项目仅作 fallback） ───────────────
# AVONET — 优先本地，否则用 ~/lucc
AVONET_PATH <- Sys.getenv(
  "V4_AVONET_PATH",
  {
    local <- file.path(DIRS$external, "traits", "AVONET1_BirdLife.csv")
    if (file.exists(local)) local else file.path("~", "lucc", "AVONET1_BirdLife.csv")
  }
)

# EltonTraits
ELTONTRAITS_PATH <- Sys.getenv(
  "V4_ELTONTRAITS_PATH",
  {
    local <- file.path(DIRS$external, "traits", "BirdFuncDat.txt")
    if (file.exists(local)) local else file.path("~", "lucc", "EltonTraits", "BirdFuncDat.txt")
  }
)

# 性状插补表（本地化优先）
TRAIT_IMPUTED_PATH <- Sys.getenv(
  "V4_TRAIT_IMPUTED_PATH",
  {
    local <- file.path(DIRS$external, "traits", "table_species_traits_imputed.csv")
    if (file.exists(local)) {
      local
    } else {
      # fallback：邻近项目
      file.path(PROJECT_ROOT, "..", "bird_grid_community_analysis", "results",
                "table_species_traits_imputed.csv")
    }
  }
)

# IUCN
IUCN_API_KEY <- Sys.getenv("IUCN_API_KEY", "")

# ── 随机森林参数 ───────────────────────────────────────────────────
RF_NUM_DRAWS    <- 100L
RF_NUM_TREES    <- 1000L
RF_SEED         <- 2025L
RF_MIN_NODE     <- 5L
RF_MTRY_FACTOR  <- NULL

# ── 驱动变量分组 ───────────────────────────────────────────────────
DRIVER_GROUPS_TREND <- list(
  climate_change   = c("delta_t_mean", "delta_t_extreme"),
  landuse_change   = c("delta_forest", "delta_impervious", "delta_natural"),
  human_change     = c("delta_hfi"),
  baseline         = c("bio11", "elev_mean", "centroid_lon", "centroid_lat")
)

DRIVER_GROUPS_TREND_SENSITIVITY <- list(
  climate_change   = c("delta_t_std"),
  landuse_change   = c("delta_forest", "delta_impervious", "delta_natural"),
  human_change     = c("delta_hfi"),
  baseline         = c("bio11", "elev_mean", "centroid_lon", "centroid_lat")
)

DRIVER_GROUPS <- list(
  climate      = c("bio4", "bio7", "bio11", "bio13"),
  topo_habitat = c("elev_mean", "elev_sd", "texture_shannon",
                   "habitat_diversity_shannon"),
  human        = c("hfi_mean", "landcover_built", "landcover_cropland"),
  space        = c("centroid_lon", "centroid_lat")
)

DRIVER_TREND_LABELS <- c(
  delta_t_mean     = "Mean temperature change",
  delta_t_extreme  = "Extreme temperature change",
  delta_t_std      = "Standardized temperature change",
  delta_forest     = "Forest cover change",
  delta_impervious = "Impervious surface change",
  delta_natural    = "Natural land cover change",
  delta_hfi        = "Human footprint change",
  bio11            = "Mean temperature (baseline)",
  elev_mean        = "Elevation (baseline)",
  centroid_lon     = "Longitude",
  centroid_lat     = "Latitude"
)

DRIVER_GROUP_TREND_LABELS <- c(
  climate_change = "Climate change",
  landuse_change = "Land use change",
  human_change   = "Human pressure change",
  baseline       = "Spatial baseline"
)

# ── 外部栅格目录 ───────────────────────────────────────────────────
CRU_DIR  <- Sys.getenv("V4_CRU_DIR",  file.path(DIRS$external, "cru_ts"))
CLCD_DIR <- Sys.getenv("V4_CLCD_DIR", file.path(DIRS$external, "clcd"))
HFI_DIR  <- Sys.getenv("V4_HFI_DIR",  file.path(DIRS$external, "hfi"))
WORLDCLIM_DIR <- Sys.getenv("V4_WORLDCLIM_DIR",
  file.path(DIRS$external, "worldclim"))

# ── Nature/Science 图表规范 ────────────────────────────────────────
NATURE_FONT      <- "Arial"
NATURE_PT        <- 7.5
NATURE_WIDTH_S   <- 89    # mm, 单栏
NATURE_WIDTH_M   <- 120   # mm, 1.5 栏
NATURE_WIDTH_L   <- 183   # mm, 双栏
NATURE_ACCENT    <- "#0E5A78"
NATURE_DPI       <- 300L
NATURE_LINE_AXIS <- 0.25
NATURE_LINE_TICK <- 0.15
NATURE_LINE_MAP  <- 0.12

# ── 地图硬规则（v4 强制） ────────────────────────────────────────────
MAP_BBOX_XLIM <- c(73, 135)
MAP_BBOX_YLIM <- c(18, 54)

# ── 检测协变量公式（v4 沿用 v3，但描述与稿件统一） ────────────────────
# log_events：log1p(visit 次数)
# log_duration：log10(mean_duration_min)，缺失时为 0（由 has_duration=0 标记）
# has_duration：1 表示有 duration 记录，0 表示缺失（missingness indicator）
DET_FORMULA_STR <- "~ log_events + log_duration + has_duration"

# ── 占有率协变量公式 ──────────────────────────────────────────────
OCC_FORMULA_STR <- paste(
  "~ bio4 + bio7 + bio11 + bio13 +",
  "elev_mean + elev_sd + texture_shannon + habitat_diversity_shannon +",
  "hfi_mean + landcover_built + landcover_cropland +",
  "centroid_lon + centroid_lat + year_scaled"
)

# ── 输出目录确保函数 ──────────────────────────────────────────────
ensure_v4_dirs <- function() {
  for (d in DIRS[c("derived", "results", "figures", "logs")]) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  DIRS
}

message("[00_config_v4] loaded  —  project: ", PROJECT_ROOT)
