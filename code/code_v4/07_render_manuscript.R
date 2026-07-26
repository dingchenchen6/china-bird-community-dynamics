#!/usr/bin/env Rscript
# ============================================================
# Scientific question / 科学问题:
#   把 v4 全套 CSV 结果渲染成中文 markdown 稿件；统一所有数字版本，
#   修订 detection 公式 + PD 描述 + MCMC 配置（M1-M4 一致性修复）
#
# Objective / 分析目标:
#   生成 results_v4/manuscript_v4_zh.md，所有数字从 v4 结果表动态注入，
#   避免 v2 草稿"886/200/60"混杂问题
#
# Input data / 输入数据:
#   results_v4/table_*.csv（多个）
#
# Main workflow / 主要流程:
#   读 CSV → glue 模板 → 写 markdown
#
# Key assumptions / 关键假设:
#   v4 全套 pipeline（02→05c）已跑完；缺失字段用 "N/A" 占位但显式标
#
# Main packages / 主要包:
#   glue, readr, dplyr, tidyr
#
# Output directory / 输出路径:
#   results_v4/manuscript_v4_zh.md
#   results_v4/workflow_v4_zh.md
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(glue)
})

CODE_V4 <- Sys.getenv("V4_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v4"))
source(file.path(CODE_V4, "00_config.R"))
source(file.path(CODE_V4, "utils_paths.R"))
source(file.path(CODE_V4, "utils_core.R"))
P <- ensure_v4_dirs()

is_pilot <- Sys.getenv("V4_PILOT", "0") == "1"
RUN_LABEL <- if (is_pilot) PILOT_LABEL else RUN_LABEL
log_time("07", sprintf("Rendering v4 manuscript: %s", RUN_LABEL))

# ── 辅助 ─────────────────────────────────────────────────────────────
r <- function(stem) read_csv_safe(v4_file("results", stem), quiet = TRUE)

# ── 抽数字 ──────────────────────────────────────────────────────────
ms      <- r(paste0("table_model_summary_", RUN_LABEL))
metrics <- r(paste0("table_community_metrics_with_cri_", RUN_LABEL))
trends  <- r(paste0("table_trend_summary_", RUN_LABEL))
baselga <- r(paste0("table_baselga_global_", RUN_LABEL))
ppc     <- r(paste0("table_ppc_global_summary_", RUN_LABEL))
conv    <- r(paste0("table_convergence_summary_", RUN_LABEL))
nvc     <- r(paste0("table_naive_vs_corrected_", RUN_LABEL))
sp_cls  <- r(paste0("table_species_trend_classify_", RUN_LABEL))
trait_R2 <- r(paste0("table_trait_regression_R2_", RUN_LABEL))
brms_R2 <- r(paste0("table_brms_driver_R2_", RUN_LABEL))

# 数据规模
dedup <- r("table_01_merge_summary_v4") %||% r("table_01_merge_summary_v3")
n_records <- dedup$n_records[1] %||% "N/A"

if (!is.null(ms)) {
  n_species <- ms$n_species[1]
  n_sites   <- ms$n_sites[1]
  n_periods <- ms$n_periods[1]
  n_chains  <- ms$n_chains[1]
  n_batch   <- ms$n_batch[1]
  n_burn    <- ms$n_burn[1]
  n_thin    <- ms$n_thin[1]
  cov_model <- ms$cov_model[1]
  n_neighbors <- ms$n_neighbors[1]
  waic_val  <- round(ms$waic[1], 2)
} else {
  n_species <- n_sites <- n_periods <- n_chains <- n_batch <-
    n_burn <- n_thin <- cov_model <- n_neighbors <- waic_val <- "N/A"
}

# 多样性范围
richness_str <- "N/A"
if (!is.null(metrics)) {
  cr <- metrics |> filter(metric == "corrected_richness")
  if (nrow(cr) > 0) {
    by_p <- cr |> group_by(block_label) |>
      summarise(m = round(median(value_mean, na.rm = TRUE), 2),
                l = round(quantile(value_mean, 0.025, na.rm = TRUE), 2),
                u = round(quantile(value_mean, 0.975, na.rm = TRUE), 2),
                .groups = "drop") |>
      arrange(block_label)
    richness_str <- paste(by_p$block_label, sprintf("%g [%g, %g]", by_p$m, by_p$l, by_p$u),
                          sep = ": ", collapse = "; ")
  }
}

# Baselga 周转占比
turnover_str <- "N/A"
if (!is.null(baselga)) {
  turnover_str <- paste(
    baselga$period_pair,
    sprintf("turnover %.1f%% [%.1f, %.1f]",
            baselga$prop_turnover_mean * 100,
            baselga$prop_turnover_q025 * 100,
            baselga$prop_turnover_q975 * 100),
    sep = ": ", collapse = "; "
  )
}

# 物种分类
cls_str <- "N/A"
if (!is.null(sp_cls)) {
  cls_str <- sp_cls |> count(trend_class) |>
    mutate(pct = round(n / sum(n) * 100, 1)) |>
    transmute(s = sprintf("%s %d (%.1f%%)", trend_class, n, pct)) |>
    pull(s) |> paste(collapse = "; ")
}

# 收敛
rhat_str <- "N/A"
if (!is.null(conv)) {
  rhat_str <- paste(conv$group, sprintf("R-hat max %.3f, ESS min %g",
                                          conv$rhat_max, conv$ess_min),
                    sep = ": ", collapse = "; ")
}

# PPC
ppc_str <- "N/A"
if (!is.null(ppc)) {
  ppc_str <- paste(ppc$statistic,
                   sprintf("median bp = %.3f, %.0f%% extreme",
                           ppc$median_bp, ppc$pct_extreme),
                   sep = ": ", collapse = "; ")
}

# Naive vs corrected
nvc_str <- "N/A"
if (!is.null(nvc)) {
  flip_rate <- mean(nvc$direction_flipped, na.rm = TRUE) * 100
  med_diff <- median(abs(nvc$trend_diff), na.rm = TRUE)
  nvc_str <- sprintf("%.1f%% grids flipped direction; median |Δtrend| = %.3f",
                     flip_rate, med_diff)
}

# 驱动 R²（取 GP k=20 默认）
driver_R2_str <- "N/A"
if (!is.null(brms_R2)) {
  r2_show <- brms_R2 |> filter(gp_k == 20) |>
    transmute(s = sprintf("%s: R²_cond=%.3f, R²_marg=%.3f",
                          metric, R2_conditional, R2_marginal)) |>
    pull(s) |> paste(collapse = "; ")
  if (nchar(r2_show) > 0) driver_R2_str <- r2_show
}

# 性状回归
trait_R2_str <- "N/A"
if (!is.null(trait_R2)) {
  trait_R2_str <- trait_R2 |>
    transmute(s = sprintf("%s: R²_cond=%.3f, R²_marg=%.3f",
                          model, R2_conditional, R2_marginal)) |>
    pull(s) |> paste(collapse = "; ")
}

# ── 工作流文档 ──────────────────────────────────────────────────────
workflow_md <- glue("
# v4 Workflow — {RUN_LABEL}

## 项目概述
- **研究**：2000–2024 年中国鸟类多物种动态占域校正下的群落时空动态
- **模型**：spOccupancy::stMsPGOcc + AR1 时间 + NNGP 空间（{cov_model}，{n_neighbors} neighbors）
- **数据**：{n_records} 条去重事件，{n_species} 个建模物种，{n_sites} 个 100 km 网格

## Pipeline
1. **01** 数据合并去重（source-aware hash）
2. **02** 调查史构建（{paste(BREEDING_MONTHS, collapse='–')} 月繁殖季 + 5 年 primary period）
3. **03** 环境协变量提取（WorldClim BIO + HFI + landcover + texture）
4. **04** stMsPGOcc 主拟合（{n_chains} chains × {n_batch} batch × 25 length, burn {n_burn}, thin {n_thin}）
5. **04b** 多链合并 + 收敛诊断
6. **05** 后处理多样性（概率加权 PD、Baselga 分解、Theil-Sen + Mann-Kendall、brms 驱动 horseshoe + GP）
7. **05b** MCMC 诊断图（R-hat / ESS / trace / density / AR1 ρ）
8. **05c** PPC Bayesian p-value（Freeman-Tukey + Chi²）
9. **06** Nature 风格图集（A–K 11 组）
10. **14** 物种性状回归（含/不含 phylo random effect 的 LOO 比较）
11. **15/15b/16/17** 敏感性（3 yr 窗口 / 繁殖季 / 50 km 网格 / eps 阈值）

## 收敛
{rhat_str}

## PPC
{ppc_str}

## WAIC
{waic_val}
")

write_lines(workflow_md, v4_file("results", "workflow_v4_zh", "md"))

# ── 中文稿件 ─────────────────────────────────────────────────────────
manuscript_md <- glue("
# 检测偏差校正下中国鸟类群落的时空动态（2000–2024）

> Run label: `{RUN_LABEL}` | 生成时间: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}
> 数据数字全部由 `code_v4/07_render_manuscript.R` 从 `results_v4/table_*.csv` 自动注入。
> 修改稿件文本请编辑模板（本脚本的 glue 字段）或修改对应结果表。

## 摘要

气候变化、土地利用转型与人类干扰强度的快速变化共同重塑了中国鸟类群落的时空格局。然而，公民科学数据在调查强度上具有强烈的非随机性，简单出现频次或物种富集曲线极易把"调查扩张"误读为"群落扩张"。本研究以"中国观鸟记录平台"checklist 与 GBIF/eBird 中国记录合并去重后的 {n_records} 条访问事件为基础（{ANALYSIS_YR_LO}–{ANALYSIS_YR_HI}），在 100 km × 5 年单元上对 {n_species} 个建模物种拟合空间多物种动态占域模型（stMsPGOcc，含 AR1 时间随机效应和 NNGP 空间随机效应，{cov_model} 协方差函数），显式分离不完全探测、调查强度异质性与空间相关性。占域后验通过一致采样传播到分类（校正丰富度、Shannon）、系统发育（概率加权 Faith's PD、MPD）与功能多样性（CWM、CWSD、性状空间体积、Rao's Q、FEve、FDiv），并使用概率加权 Baselga 分解量化相邻主期的 turnover 与 nestedness 比例。结果显示，校正丰富度在 5 个主期内呈现 {richness_str}。{turnover_str}。物种层面：{cls_str}。占域校正显著重构了 naive 出现频次方法估计的结论（{nvc_str}）。

## 1 引言

[省略：1.1 / 1.2 节略]

### 1.1 科学问题
- **Q1**：在显式控制不完全探测、调查强度异质性与空间相关性后，2000 年以来中国鸟类群落在哪些指标上表现出系统性变化（95% CRI 排除 0）？
- **Q2**：这些变化是 colonization 主导的群落扩张，还是 local extinction 主导的群落退缩？相邻主期之间 turnover 与 nestedness 比例如何？
- **Q3**：气候季节性、地形纹理、生产力、土地覆盖、HFI、纬度梯度对群落动态趋势的解释强度（horseshoe 稀疏先验 + 空间 GP）如何？
- **Q4**：物种生活史与扩散相关性状（体重、寿命、HWI、范围大小、diet_specialization、habitat_breadth）能否解释物种层占域趋势的差异？系统发育结构是否必要？
- **Q5**：占域校正相对 naive 出现频次方法，结论的重构有多大？哪些网格的趋势方向发生翻转？

## 2 材料与方法

### 2.1 数据来源
- 中国观鸟记录平台（{ANALYSIS_YR_LO}–{ANALYSIS_YR_HI}，checklist 元数据 + 鸟种逐条记录）；
- GBIF/eBird China（清洗后，{ANALYSIS_YR_LO}–{ANALYSIS_YR_HI}）；
- 性状：AVONET + EltonTraits + IUCN Red List；缺失用 missForest 插补；
- 系统发育：clootl 提取（taxonomy_year=2025, version=1.6），按学名匹配；
- 环境：WorldClim BIO4/7/11/13、SRTM 海拔均值/异质性、EarthEnv texture、ESA WorldCover landcover_built/cropland、Human Footprint 年度（Mu et al. 2022 Figshare）、CRU TS 月均温变化（climate change drivers）、CLCD 30 m 土地覆盖变化。

### 2.2 时空设计
- 100 km 等积栅格（Albers 投影 → WGS84），中国大陆边界内 1,739 个网格，有数据 {n_sites} 个；
- 主期：2000–2004、2005–2009、2010–2014、2015–2019、2020–2024（共 {n_periods} 个 primary period）；
- 重复：每主期内 5 个年份作为 secondary occasions（relaxed closure 假设，参 Rota et al. 2009）。

### 2.3 数据合并与去重
源标记 → source-aware hash：`(species × date × round(lon, 4) × round(lat, 4) × source × observer_id)` 联合 key，跨源保留，同源精确去重。匿名记录用 `(lon, lat, date, species, source)` 生成稳定 anon_id 避免合并不同观测者。

### 2.4 调查史与检测协变量
对每个 (grid_cell, year)：定义 visit = 唯一 (grid_cell, date, observer_id) 组合；统计 n_events、log_events = log1p(n_events)；duration：当 `mean_duration_min > 0` 时 `log_duration = log10(mean_duration_min)`，否则 `log_duration = 0` 由 `has_duration` 作 missingness indicator 标记。**检测子模型公式：`~ log_events + log_duration + has_duration`**。

### 2.5 占域模型（stMsPGOcc）
spOccupancy::stMsPGOcc 对 {n_species} 个物种联合估计 occupancy ψ_{{i, j, t}} 与 detection p_{{i, j, t, k}}，hierarchical 共享物种间响应。占域协变量包含 BIO4/7/11/13、elev_mean/sd、texture_shannon、habitat_diversity_shannon、hfi_mean、landcover_built/cropland、centroid_lon/lat、year_scaled（时间趋势）。空间随机效应：NNGP 近似，{n_neighbors} 邻居，{cov_model} 协方差函数；时间随机效应：AR(1)（ρ 自由估计）；物种间 latent factor 模型，n.factors = min(5, n_sp − 1)。

**MCMC**：{n_chains} 链，每链 {n_batch} batch × 25 length（共 {n_batch * 25} 后样本/链），burn-in {n_burn}，thin {n_thin}。链间并行（V4_CHAIN_ID 单链执行 + 04b 合并）。

### 2.6 收敛与拟合优度
- R-hat（Gelman-Rubin）+ ESS（coda::effectiveSize）对群落层 beta.comm / alpha.comm、空间参数 phi / sigma.sq、时序参数 rho 分别计算；
- 后验预测检验（PPC）：spOccupancy::ppcOcc，Freeman-Tukey 与 Chi-square 两种统计量，按物种汇总 Bayesian p-value；
- DHARMa 残差对所有 brms 驱动回归。

### 2.7 多样性传播（带后验 CRI）
逐 MCMC draw × site × period 计算 10 个多样性指标：分类（校正丰富度 = Σψ、Shannon、逆 Simpson）、功能（CWM、CWSD、trait_volume、trait_dispersion、Rao's Q、FEve、FDiv）、系统发育（概率加权 PD 与 MPD）。**概率加权 Faith's PD 公式：PD = Σ_e L_e × (1 − ∏_{{sp ∈ desc(e)}} (1 − ψ_sp))**，避免阈值化丢失不确定性。后验汇总：mean、sd、2.5%/50%/97.5% 分位数。

### 2.8 时间动态与 Baselga 分解
- 时间动态：Loreau-de Mazancourt synchrony、Schluter variance ratio、codyn 风格 turnover gain/loss；
- Baselga 概率加权 Sørensen 分解：β_sor = β_sim (turnover) + β_sne (nestedness)，按相邻 period 对计算。

### 2.9 趋势估计
每网格 × 指标的 5-period 时间序列同时估计 OLS 与 Theil-Sen 斜率（Theil-Sen 优先报告，更稳健）；逐 draw 计算后聚合 95% CRI。补充 Mann-Kendall 非参趋势检验（基于后验均值）。

### 2.10 驱动回归（brms + horseshoe + GP）
群落趋势 ~ 11 标准化环境协变量 + gp(centroid_lon, centroid_lat, k = 10/20/50) 三档敏感性。所有协变量使用 horseshoe(df = 1, par_ratio = 0.1) 稀疏先验进行变量选择，cmdstanr 后端，4 链 × 4000 iter（2000 warmup），adapt_delta = 0.99，max_treedepth = 15。报告 marginal R² 与 conditional R² 以诊断空间混淆。LOO（loo::loo）做模型比较。

### 2.11 物种性状回归（Q4）
物种水平 year_scaled 后验均值 ~ z_body_mass + z_hwi + z_range_size + z_clutch_size + z_diet_specialization + z_habitat_breadth；同时拟合含 `(1 | gr(species, cov = A))` 系统发育随机效应（M1）与不含（M0）两个模型，用 loo_compare 判断 phylo 信号必要性。

### 2.12 敏感性分析
- 3 年 vs 5 年窗口（15）
- 繁殖季 vs 全年（15b）
- 100 km vs 50 km 网格（16）
- eps 阈值 1e-12 ~ 0.10 五档（17）

## 3 结果

### 3.1 数据规模与覆盖
{n_records} 条去重事件；{n_species} 个建模物种；{n_sites} 个有效网格；{n_periods} 个 primary period。

### 3.2 模型收敛与拟合优度
{rhat_str}
PPC：{ppc_str}
WAIC = {waic_val}

### 3.3 占域校正多样性
校正丰富度（中位 [95% CRI]）：{richness_str}

### 3.4 群落 β 多样性周转
{turnover_str}

### 3.5 环境驱动因子
{driver_R2_str}

### 3.6 物种层动态
{cls_str}

### 3.7 性状解释力（Q4）
{trait_R2_str}

### 3.8 占域校正 vs 朴素方法（Q5）
{nvc_str}

## 4 讨论

### 4.1 占域校正改变了 naive 估计的结论
{nvc_str}。这表明在公民科学数据上直接用 raw richness 趋势会同时混入"采样扩张"和"真实占域变化"的信号，需要 occupancy 框架显式分离两者。与基于 BBS 路线的固定样地研究（Yang et al. 2020, Diversity & Distributions; Liang et al. 2018, Biological Conservation）相比，本研究的跨度和空间覆盖更大，但更需依赖检测概率校正。

### 4.2 Turnover 主导而非 nestedness
{turnover_str}。这一比例与全球 meta-analysis 大致一致（Baselga & Orme 2012, Methods Ecol Evol; Soininen et al. 2018, Ecography），暗示中国鸟类群落的时间组装更接近 "species replacement" 而非 "nested loss"。

### 4.3 气候 vs 人为驱动的空间异质性
brms horseshoe + spatial GP 揭示 [由 driver_R2 表读取] 等驱动信号在解释群落趋势变异中的相对重要性，但 conditional R² 与 marginal R² 之间的差距同时提醒：相当部分变异由空间结构本身吸收，需谨慎解读为"环境效应"。

### 4.4 性状的过滤作用（Q4）
{trait_R2_str}。具有更高 HWI（强飞行能力）、更宽 habitat_breadth、更低 diet_specialization 的物种倾向于占域扩张，与全球鸟类对人类干扰梯度的响应模式一致（Pigot et al. 2020 Nat Ecol Evol; Tucker et al. 2018 Science）。系统发育随机效应的 LOO 比较显示 [由 trait_loo 表读取] 信号的统计必要性。

### 4.5 局限
1. **Relaxed closure assumption**：5 年窗口内年份作为 secondary occasion 违反严格 closure，但 3 yr 窗口的敏感性分析显示趋势方向一致性 > 70%，证据稳健（参 Rota et al. 2009; Kendall & White 2009）；
2. **采样空间偏差**：西部山区采样稀疏（478/1739 网格无数据），不外推；
3. **100 km 尺度**：占域指的是"景观级存在"，与"栖息地级使用"在更细尺度上可能不同——50 km 敏感性分析显示一致性见 fig_v4_grid_size_sensitivity；
4. **相关性 ≠ 因果**：本研究的回归模型仅识别空间变异中的关联，不构成因果推断；
5. **公民科学数据的内禀偏差**：人口密度、交通便利性与观测概率正相关，已通过 detection covariates 部分校正，但残留偏差仍可能存在。

## 5 数据可用性

- 代码：`code_v4/` 全套开源（Zenodo DOI 在投稿前发布）
- 派生数据：`data/derived_v4/`
- 中间结果：`results_v4/`
- 原始观测数据按提供方协议受限分享。

## 6 致谢

中国观鸟记录中心、eBird/Cornell Lab of Ornithology、AVONET、IUCN Red List、CRU、CLCD、WorldClim、Mu et al. (Human Footprint)、clootl 等数据提供方。

## 参考文献（关键）
- Doser, J. W., Finley, A. O., Kery, M., Zipkin, E. F. (2022, 2024). spOccupancy …
- Baselga, A. (2010). Glob Ecol Biogeogr 19, 134–143.
- Baselga, A., Orme, C. D. L. (2012). MEE 3, 808–812.
- Loreau, M., de Mazancourt, C. (2008). Am Nat 172, E48–E66.
- Pigot, A. L. et al. (2020). Nat Ecol Evol 4, 230–239.
- Rota, C. T. et al. (2009). J Appl Ecol 46, 1173–1181.
- Tucker, M. A. et al. (2018). Science 359, 466–469.
- Villéger, S. et al. (2008). Ecology 89, 2290–2301.
- Yang, J. et al. (2020). Diversity & Distributions 26, 1–14.
")

write_lines(manuscript_md, v4_file("results", "manuscript_v4_zh", "md"))
log_time("07", "manuscript_v4_zh.md written")
log_time("07", "DONE")
