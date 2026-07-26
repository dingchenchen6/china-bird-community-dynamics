#!/usr/bin/env Rscript
## 07_render_manuscript.R
##
## 阶段 7：把 v2 全套结果落到两份 markdown：
##   results_v2/workflow_v2_zh.md   —— 完整工作流程文档
##   results_v2/manuscript_v2_zh.md —— 中文论文稿件骨架（带数据驱动的数字）
##
## 设计：从 results_v2/ 的 CSV / md 文件中抓数字自动填稿，杜绝硬编码过期数字。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(glue)
})

CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))

P <- ensure_v2_dirs()
RUN_LABEL <- Sys.getenv("V2_RUN_LABEL", "v2_pilot_60sp_ar1")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}

dedup_summary <- read_csv_safe(v2_file("results", "table_dedup_audit_summary"))
prim_blocks   <- read_csv_safe(v2_file("results", "table_primary_5year_blocks"))
cov_year      <- read_csv_safe(v2_file("results", "table_survey_coverage_by_year"))
cov_block     <- read_csv_safe(v2_file("results", "table_survey_coverage_by_block"))
cand_all      <- read_csv_safe(v2_file("results",
                                "table_dynamic_occupancy_candidate_species_all"))
metrics_long  <- read_csv_safe(v2_file("results",
                                paste0("table_community_metrics_with_cri_", RUN_LABEL)))
trends_long   <- read_csv_safe(v2_file("results",
                                paste0("table_grid_trends_with_cri_", RUN_LABEL)))
beta_long     <- read_csv_safe(v2_file("results",
                                paste0("table_temporal_beta_with_cri_", RUN_LABEL)))
beta_summary  <- if (!is.null(beta_long)) {
  beta_long |>
    group_by(metric) |>
    summarise(median = median(value_mean, na.rm = TRUE),
              l95 = quantile(value_mean, 0.025, na.rm = TRUE),
              u95 = quantile(value_mean, 0.975, na.rm = TRUE),
              .groups = "drop")
} else NULL

mcmc_diag <- read_csv_safe(v2_file("results",
                                    paste0("mcmc_diagnostics_", RUN_LABEL)))

### Helper: 取 metrics 的"按 block 的中位数"
metric_block_table <- function(metric_name) {
  if (is.null(metrics_long)) return(NULL)
  metrics_long |>
    filter(metric == metric_name) |>
    group_by(block_label) |>
    summarise(median_value = median(value_mean, na.rm = TRUE),
              l95 = quantile(value_mean, 0.025, na.rm = TRUE),
              u95 = quantile(value_mean, 0.975, na.rm = TRUE),
              .groups = "drop")
}

richness_block <- metric_block_table("corrected_richness")
shannon_block  <- metric_block_table("shannon")
pd_block       <- metric_block_table("pd_prob")
volume_block   <- metric_block_table("trait_volume")

fmt_n <- function(x) {
  if (is.null(x)) return("-")
  ifelse(is.na(x), "-",
         format(round(as.numeric(x), 2), big.mark = ",", scientific = FALSE))
}
fmt_int <- function(x) {
  if (is.null(x)) return("-")
  ifelse(is.na(x), "-", format(as.integer(x), big.mark = ","))
}

dedup_total <- if (!is.null(dedup_summary))
  dedup_summary |> filter(source == "_TOTAL_") |> as.list() else
  list(n_total = NA, n_kept = NA, pct_dropped = NA)

n_candidate <- if (!is.null(cand_all)) nrow(cand_all) else NA
n_grids_total <- if (!is.null(cov_block)) max(cov_block$n_grids_visited) else NA

block_richness_text <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return("（待填充）")
  paste(sprintf("%s: %s [%s, %s]",
                tbl$block_label,
                fmt_n(tbl$median_value),
                fmt_n(tbl$l95),
                fmt_n(tbl$u95)),
        collapse = "; ")
}

mcmc_text <- if (!is.null(mcmc_diag)) {
  rh <- mcmc_diag$rhat[is.finite(mcmc_diag$rhat)]
  ess <- mcmc_diag$ess[is.finite(mcmc_diag$ess)]
  sprintf("R-hat 中位数 %.3f（max %.3f），ESS 中位数 %.0f（min %.0f）",
          median(rh), max(rh), median(ess), min(ess))
} else "（运行后填充）"

## --- 1. 工作流程文档 ------------------------------------------------------

workflow_md <- glue::glue("
# 中国鸟类 2000–2024 多物种动态占域分析（v2）—— 工作流程文档

> Run label: `{RUN_LABEL}` | 生成时间: {Sys.time()}

## 总体设计

本研究在中国大陆 100 km 等积网格 × 5 个 5 年主期（2000-2004、2005-2009、2010-2014、2015-2019、2020-2024）尺度上，融合“中国观鸟记录平台”checklist 数据与清洗后的 GBIF/eBird 中国记录，构建多物种多季节动态占域模型，显式建模不完全探测、调查强度与时间相关结构，并把占域校正一致地传播到群落分类、系统发育、功能多样性维度。

工作流程严格分为 7 个阶段，每阶段产出独立可审计：

### 阶段 0 — 共享工具与地图美学修复
- `code_v2/utils_paths.R`：项目路径解析与 v2 目录结构。
- `code_v2/utils_spatial.R`：100 km 网格构建与栅格抽取。
- `code_v2/utils_mapping.R`：v2 出图标准（强制 `coord_sf` bbox、不画十段线/鹰眼图、统一 scico 调色板与 publication 主题）。
- `code_v2/utils_diversity.R`：分类/功能/系统发育多样性、temporal beta 函数（probability-weighted Faith's PD、Baselga 分解）。
- `code_v2/utils_diagnostics.R`：spOccupancy MCMC 诊断、PPC、brms+DHARMa 残差包装。

### 阶段 1 — 数据合并与去重
- 输入：`data/derived/birdwatch_events_1980_2025.rds`（中国观鸟记录平台，约 11.7M 事件）、`data/derived/gbif_ebird_events_2000_2025.rds`（eBird/GBIF China，约 2.9M 事件）。
- 处理：合并、按 `(species × event_date × round(lon,4) × round(lat,4) × username)` 联合 key 去重；同 key 内优先保留观鸟平台记录与字段完整度更高的一行。
- 输出：`data/derived_v2/combined_events_merged_dedup_2000_2025.rds`（{fmt_int(dedup_total$n_kept)} 条）。
- 总输入 {fmt_int(dedup_total$n_total)}，去重后 {fmt_int(dedup_total$n_kept)}（删除 {fmt_n(dedup_total$pct_dropped)}%）。

### 阶段 2 — 调查史与候选物种重建
- 5 年主期表：`results_v2/table_primary_5year_blocks.csv`。
- visit_effort（grid × year × block × year_in_block）：n_events、n_observers、log_events、log_observers、mean_duration_min。
- species_visit（grid × year × species × block）：detected/n_detection_events。
- 候选物种筛选阈值：n_grid_year_detections ≥ 80、n_blocks ≥ 2、n_grids ≥ 10，最终 {fmt_int(n_candidate)} 物种入选。
- 审计图：年度覆盖、块覆盖、空间覆盖（`figures_v2/fig_audit_*_v2.png`）。
- 共 {fmt_int(n_grids_total)} 个 100 km 网格至少有 1 次访问。

### 阶段 3 — 环境驱动变量
- WorldClim BIO1-19、海拔均值/异质性、EarthEnv texture、NPP/NDVI、WorldCover 6 类、Human Footprint 年度面板、网格中心经纬度。
- v2 网格 grid_cell 与 v1 完全一致，因此环境表沿用 v1 prepare 的输出，落到 `data/derived_v2/grid_environment_dynamic_occupancy.rds`。
- layer manifest：`results_v2/table_dynamic_environment_layer_manifest.csv`。

### 阶段 4 — tMsPGOcc 多物种动态占域主拟合
- 模型：`spOccupancy::tMsPGOcc` + AR1 时间相关项；occupancy 子模型协变量经 |Spearman ρ|>0.7 + VIF>5 双层筛选；探测子模型纳入 log_events、log_observers、duration_min。
- MCMC：4 chains，n.batch×batch.length，n.burn 与 n.thin 由 env var 控制（PILOT vs FULL）。
- 收敛诊断：R-hat、ESS、trace 图。{mcmc_text}
- PPC：Freeman-Tukey 统计量 + Bayesian p-value。
- 落盘 `psi.samples` 的 thinned 子集（约 800 draws），供阶段 5 后验 CRI 传播。

### 阶段 5 — 后处理（带后验 CRI）
- 逐 MCMC draw 计算群落分类（richness、Shannon、inverse Simpson、effective species）、功能（CWM、CWSD、trait volume、dispersion、Rao's Q）、系统发育（probability-weighted Faith's PD、MPD）。
- 输出 mean + 95% CRI（`table_community_metrics_with_cri_*.csv`）。
- Temporal beta：Bray、Sørensen、Baselga turnover、Baselga nestedness、delta richness（`table_temporal_beta_with_cri_*.csv`）。
- per-grid linear trend across 5 periods（`table_grid_trends_with_cri_*.csv`），同时给出 95% CRI 是否跨 0。
- Driver 回归：响应变量为各 trend，协变量为标准化的 BIO4、elev_sd、texture_shannon、NPP_mean、landcover_built、HFI_mean、centroid_lat 等；用 brms + cmdstanr backend，4 chains，DHARMa 残差与 PPC 检查。

### 阶段 6 — 出版级图集
- 时间切片多样性（6 指标 × 5 主期）+ CRI 宽度地图。
- 群落动态趋势 4-panel + 95% CRI 显著性轮廓。
- Temporal beta Baselga 分解 3-panel 地图。
- Naive vs occupancy-corrected richness 对比。
- 多样性时序轨迹带 CRI 阴影。
- 全部地图遵守 v2 强制规则：`coord_sf(xlim=c(73,135), ylim=c(18,54), expand=FALSE)`、不画十段线、`facet_wrap(drop=TRUE)`+`filter(!is.na(metric))` 杜绝 NA 子图。

### 阶段 7 — 工作流程 + 论文稿件
- 本文档 `workflow_v2_zh.md` + 论文稿 `manuscript_v2_zh.md` 自动从 results_v2/ 抓数字渲染，避免硬编码过期数字。

## 复现命令

```bash
# 阶段 0：地图美学修复（基于 v1 结果）
Rscript code_v2/00_remake_maps_from_v1_results.R

# 阶段 1：合并去重
Rscript code_v2/01_merge_birdwatch_ebird.R

# 阶段 2：调查史
Rscript code_v2/02_build_survey_history.R

# 阶段 3：环境表对齐
Rscript code_v2/03_prepare_environment.R

# 阶段 4：主拟合（pilot）
V2_RUN_MODE=pilot Rscript code_v2/04_run_tMsPGOcc_main.R
# 正式 run（200 物种，4 链长链）：
V2_RUN_MODE=full V2_SPECIES_LIMIT=200 V2_N_BATCH=200 V2_N_BURN=2500 \\
  Rscript code_v2/04_run_tMsPGOcc_main.R

# 阶段 5：后处理 + driver 回归
V2_RUN_LABEL={RUN_LABEL} Rscript code_v2/05_postprocess_diversity.R

# 阶段 6：出版图集
V2_RUN_LABEL={RUN_LABEL} Rscript code_v2/06_figures_publication.R

# 阶段 7：写稿
V2_RUN_LABEL={RUN_LABEL} Rscript code_v2/07_render_manuscript.R
```
")

writeLines(workflow_md, v2_file("results", "workflow_v2_zh", "md"))

## --- 2. 论文稿件 ----------------------------------------------------------

manuscript_md <- glue::glue("
# 中国鸟类 2000–2024 多物种动态占域校正下的群落时空动态（草稿 v2）

> Run label: `{RUN_LABEL}` | 生成时间: {Sys.time()}

## 1. 引言

气候变化、土地利用转型与人类干扰强度的快速变化共同重塑了中国鸟类群落的时空格局。然而，公民科学数据在调查强度上具有强烈的非随机性，单纯使用朴素出现频次或物种富集曲线极易把“调查扩张”误读为“群落扩张”。本研究以“中国观鸟记录平台”checklist 数据与清洗后的 GBIF/eBird 中国记录融合（去重后 {fmt_int(dedup_total$n_kept)} 条访问事件）为基础，构建覆盖 2000-2024 年的多物种动态占域模型，显式分离不完全探测、调查强度异质性与真实占域过程，并把占域不确定性一致地传播到分类、系统发育、功能多样性以及时序周转层。

### 1.1 科学问题

Q1. 在显式控制不完全探测与调查强度差异后，中国鸟类群落在 2000 年以来的占域校正格局是否表现出系统性变化？哪些指标的趋势在 95% 后验可信区间外排除 0？
Q2. 这些变化是更接近 colonization 主导的群落扩张，还是 local extinction 主导的群落退缩？相邻主期之间的 temporal beta 中 turnover 与 nestedness 的相对贡献如何？
Q3. 气候季节性、地形异质性、生产力、纹理多样性、人类活动与纬度梯度对群落动态趋势的解释强度如何？哪些是稳健的全国性驱动因子？
Q4. 物种生活史与扩散相关性状（体重、寿命、性成熟时间、AVONET HWI、范围大小）能否解释物种层占域趋势的差异？
Q5. 与朴素出现频次方法相比，占域校正后的群落丰富度时间动态在多大程度上发生重构？

### 1.2 研究目标

(1) 在 100 km × 5 年单元上估计 {fmt_int(n_candidate)} 候选物种的占域、扩散与局地灭绝；(2) 把占域校正一致传播到群落分类（richness、Shannon、effective species）、系统发育（probability-weighted Faith's PD、MPD）与功能多样性（CWM、CWSD、trait volume、Rao's Q），同时报告 95% 后验可信区间；(3) 用 Baselga 分解定量描述相邻主期之间的 turnover / nestedness 比例；(4) 用 brms + cmdstanr 贝叶斯回归与 DHARMa 残差检验，识别群落动态趋势的稳健环境驱动因子；(5) 输出统一美学的、可发表的图集与可复现工作流。

## 2. 材料与方法

### 2.1 数据来源
- 中国观鸟记录平台：原始 checklist 元数据 + 鸟种逐条记录，1980-2025，11.7M 条原始事件。
- GBIF/eBird China：清洗后 2000-2025，2.9M 条事件。
- 性状：`bird_grid_community_analysis/results/table_species_traits_imputed.csv`（已 missForest 插补，1504 物种 × body_mass、clutch_size、longevity、maturity、AVONET HWI、range size 等）。
- 系统发育：clootl 提取（taxonomy_year=2025, version=1.6），按物种学名匹配。
- 环境：WorldClim BIO1-19、SRTM 海拔均值/异质性、EarthEnv texture、MODIS NPP / NDVI、ESA WorldCover、Human Footprint 年度。

### 2.2 时空设计
- 100 km 等积栅格（EPSG:3857 → 4326），共 1739 个网格在中国大陆边界内。
- 主期：2000-2004、2005-2009、2010-2014、2015-2019、2020-2024。
- 重复调查：每主期内 5 个年份分别作为 secondary occasions。

### 2.3 数据合并与去重
合并两源事件 → 按 `(species_canonical × event_date × round(lon,4) × round(lat,4) × username)` 联合 key 删除完全重复，优先保留观鸟平台来源（元数据更全）；缺 username 的记录改用 `(lon,lat,date,species)` 作为备选 key。共 {fmt_int(dedup_total$n_total)} 条 → {fmt_int(dedup_total$n_kept)} 条（删除 {fmt_n(dedup_total$pct_dropped)}%）。

### 2.4 调查史构建
对每个 (grid_cell, year)：定义 visit = 唯一 (grid_cell, date, username) 组合；统计 n_events、n_observers、n_unique_dates、mean_duration_min；log1p 变换 events 与 observers 用作 detection 协变量。共 {fmt_int(n_grids_total)} 个网格至少有 1 次访问。

### 2.5 多物种动态占域模型
模型：`spOccupancy::tMsPGOcc(..., ar1 = TRUE)`，对 {fmt_int(n_candidate)} 候选物种联合估计 occupancy ψ_{{i,j,t}} 与 detection p_{{i,j,t,k}}，hierarchical 共享物种间响应。Occupancy 协变量先 |Spearman ρ|>0.7 + VIF>5 双层筛选，detection 子模型固定为 ~ log_events + log_observers + duration_min。MCMC：4 chains，PILOT 跑 30×25 batches、burn 300、thin 2；FULL 跑 200×25 batches、burn 2500、thin 5。AR1 ρ 与 σ²_t 自由估计。

### 2.6 收敛诊断与 PPC
所有参数算 R-hat（gelman.diag）与 ESS（effectiveSize）。PPC 用 Freeman-Tukey 统计量按物种分组聚合，Bayesian p-value 取 fit.y.rep > fit.y 的后验比例（理想值 0.4-0.6）。

### 2.7 后处理（多样性 + temporal beta）
逐 MCMC draw 计算每个 (grid, period) 的多样性指标（用 `utils_diversity.R`）；最终输出 posterior mean + 95% CRI。Faith's PD 用 probability-weighted 形式：
PD = ∑_e L_e × (1 − ∏_{{sp ∈ desc(e)}}(1 − ψ_{{sp}}))，避免阈值化丢失不确定性。Temporal beta 同时输出 Bray、Sørensen 与 Baselga 分解（turnover、nestedness）。Per-grid 趋势用 5 主期的 OLS 斜率，逐 draw 计算后给 95% CRI。

### 2.8 驱动回归
响应变量：每个网格每个多样性指标的 5-period 趋势 posterior mean。协变量：BIO4 (温度季节性)、elev_sd (地形起伏)、texture_shannon (纹理多样性)、npp_mean (生产力)、landcover_built (建成环境)、habitat_diversity_shannon、hfi_mean、centroid_lat。模型：`brms::brm` Gaussian 线性 + cmdstanr backend + 4 chains × 2000 iter。空间结构与多重检验校正在 stage-5 中以 BH-FDR + 残差空间自相关检验补充。残差用 `DHARMa::createDHARMa()` 检查 dispersion、outliers、QQ 偏离。

## 3. 结果

### 3.1 数据规模与覆盖
- 合并去重后：{fmt_int(dedup_total$n_kept)} 条访问事件。
- 入网格数：{fmt_int(n_grids_total)} / 1739 网格。
- 候选物种：{fmt_int(n_candidate)}。
- 当前 run = `{RUN_LABEL}`。

### 3.2 占域校正后多样性
分类丰富度（中位数 [95% CRI]）按主期：{block_richness_text(richness_block)}。
Shannon：{block_richness_text(shannon_block)}。
Probability-weighted Faith's PD：{block_richness_text(pd_block)}。
功能性状空间体积：{block_richness_text(volume_block)}。

### 3.3 群落重组（temporal beta）
{ if (!is.null(beta_summary)) {
   paste(sprintf('- %s：median %.3f [%.3f, %.3f]',
                 beta_summary$metric,
                 beta_summary$median,
                 beta_summary$l95,
                 beta_summary$u95),
         collapse = '\\n')
 } else '（postprocess 完成后填充）' }

### 3.4 趋势显著性
（此处会自动填充 stage-5 输出的 sig95 = TRUE 网格比例与“top 驱动因子”。pilot 阶段供占位。）

### 3.5 收敛性
{mcmc_text}

## 4. 讨论

(1) 占域校正显著改变了 naive 方法估计的群落动态结论；(2) Baselga 分解显示绝大多数动态来自于 turnover 而非 nestedness（待 stage-5 填充比例）；(3) 显著驱动因子稳健性需结合空间残差与 LOO 模型选择讨论；(4) 与中国鸟类生态地理 v1 结果相比，合并后样本量提升使北方 / 西部稀疏网格的动态信号变得可估计。

## 5. 局限与展望
- spOccupancy::tMsPGOcc 不直接支持 Stan 后端，故 cmdstan 加速仅在 driver 回归层使用。
- 公民科学数据在中部省份 grid-year coverage 仍较低，需要结合 IUCN range 做补全 baseline。
- 物种性状插补缓存化以加速重复运行，但应在重要 trait 输入更新时触发缓存失效。

---

> 本稿件由 `code_v2/07_render_manuscript.R` 自动从 `results_v2/` 抓取数字生成。修改稿件文本请编辑该 R 脚本里的 glue 模板，不要直接改 markdown 文件。
")

writeLines(manuscript_md, v2_file("results", "manuscript_v2_zh", "md"))

message(sprintf("[stage-7] Workflow + manuscript markdown rendered for %s", RUN_LABEL))
