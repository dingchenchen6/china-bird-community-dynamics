# Supplementary Information — Outline

> 与 main text 同步组织；所有数字字段在 FULL 跑完后由 `07_render_manuscript.R` 自动注入。

---

## S1 Data flow & deduplication

### S1.1 Source merging
- Table S1.1：原始 China Birdwatching Records + eBird China 记录数（按年）
- Fig S1.1：合并/去重前后年份事件数曲线（按 source）

### S1.2 Source-aware deduplication
- Table S1.2：去重前 / 跨源保留 / 同源精确去重的删减比例
- 与 v2 比较：v2 vs v4 记录数 diff 表
- 验证：随机抽样人工核验（n = 200 对）的真重复 / 假重复比例

### S1.3 Breeding-season filtering audit
- Table S1.3：按 period × source 的繁殖季过滤前 / 后事件数与保留率
- ⚠️ 若任一 period 保留率 <50%，加 Discussion 注释
- Fig S1.3：年份 × 月份的事件密度热力图（识别越冬偏差）

---

## S2 MCMC convergence diagnostics

### S2.1 R-hat & ESS
- Fig S2.1：R-hat 直方（beta.comm / alpha.comm / phi / sigma.sq / rho 五组）+ 1.05 阈值
- Fig S2.2：ESS 直方 + 200 阈值
- Table S2.1：max R-hat / min ESS 按参数组汇总

### S2.2 Trace & posterior density
- Fig S2.3：核心 12 个参数的 4 链 trace（每参数一面板）
- Fig S2.4：同样参数的 posterior density 叠合

### S2.3 Spatial & temporal parameters
- Fig S2.5：phi（空间衰减）后验密度 → 推断空间相关尺度
- Fig S2.6：rho（AR1 时间相关）按物种密度叠合，高亮强 AR1 物种

---

## S3 Posterior predictive checks

### S3.1 Bayesian p-value
- Table S3.1：Freeman-Tukey 与 Chi-square 两种统计量的 Bayesian p-value 按物种分布（mean、SD、% extreme < 0.05 或 > 0.95）
- Fig S3.1：Bayesian p-value 直方（两种统计量）

### S3.2 Fit.y vs fit.y.rep
- Fig S3.2：top-12 探测物种的 fit.y vs fit.y.rep 散点（4 列 × 3 行）

---

## S4 Sensitivity analyses

### S4.1 3-year vs 5-year window (script 15)
- Table S4.1：每网格 corrected_richness trend 在 3yr vs 5yr 窗口下的方向一致性
- Fig S4.1：配对散点（颜色标方向翻转）

### S4.2 Breeding-season vs year-round (script 15b)
- Table S4.2：两种过滤策略下趋势方向一致性
- Fig S4.2：配对散点

### S4.3 100-km vs 50-km grid (script 16)
- Table S4.3：50 km 子网格中位 trend 与 100 km 父网格 trend 的相关
- Fig S4.3：配对散点（颜色标方向翻转）

### S4.4 eps threshold (script 17)
- Table S4.4：5 档 eps 下 trait_volume / rao_q / FEve / FDiv 的两两 Spearman ρ
- Fig S4.4：相关性热力图

### S4.5 brms GP basis size k (10/20/50)
- Table S4.5：三档 k 下 marginal vs conditional R²
- Fig S4.5：系数 forest plot 按 k 分色

### S4.6 Candidate species threshold
- Table S4.6：80/100/150 三种 detection threshold 下物种集合差异及趋势一致性

---

## S5 Full figure atlas (Nature/Science style)

按 standard_occupancy_figures.md 的 A–K 11 组，每张图列：编号、文件名（含 .pdf 和 .png）、stage、简要说明、是否在 main text。

| # | File stem | Stage | Caption | In main? |
|---|---|---|---|---|
| A1 | fig_v4_mcmc_rhat_histogram | 05b | R-hat 直方 + 1.05 阈值 | SI |
| A2 | fig_v4_mcmc_ess_histogram | 05b | ESS 直方 + 200 阈值 | SI |
| A3 | fig_v4_mcmc_trace_core_params | 05b | 12 核心参数 4 链 trace | SI |
| A4 | fig_v4_mcmc_density_core_params | 05b | 同上 posterior density 叠合 | SI |
| B1 | fig_v4_ppc_fit_y_vs_rep | 05c | top-12 物种 fit.y vs rep 散点 | SI |
| C1 | fig_v4_beta_comm_caterpillar | 06 (TBD) | 群落 beta.comm caterpillar | SI |
| D1 | fig_v4_species_covariate_heatmap | 06 | 物种 × 协变量 系数 heatmap | **Main Fig. 3 candidate** |
| E1 | fig_v4_multidiversity_timeslices | 06 | 8 metric × 5 period 时间切片 | **Main Fig. 1** |
| E2 | fig_v4_multidiversity_uncertainty | 06 | 同上但 fill = 95% CRI width | SI |
| F1 | fig_v4_mcmc_ar1_rho_density_by_species | 05b | 物种 AR1 ρ 后验密度 | SI |
| G1 | fig_v4_multidiversity_trajectories | 06 | 9 metric 时序 + CRI 阴影 | SI |
| G2 | fig_v4_community_trends_with_cri | 06 | 4 panel 趋势地图（含 95% CRI mask） | **Main Fig. 2** |
| H1 | fig_v4_temporal_beta_baselga | 06 | Baselga 双面板地图 | SI |
| H2 | fig_v4_baselga_proportion | 06 | 全局周转占比柱状 | **Main Fig. 4** |
| I1 | fig_v4_brms_driver_forest | 06 | brms 系数 forest（含 GP k 三档） | **Main Fig. 5** |
| I2 | dharma_driver_* | 05 | DHARMa 4 panel | SI |
| J1 | fig_v4_trait_regression_forest | 14 (TBD) | 性状回归 forest | SI |
| K1 | fig_v4_naive_vs_corrected_trend | 06 | 朴素 vs 校正配对图 | **Main Fig. 6** |
| K2 | fig_v4_species_hotspot_map | 06 | 扩张/收缩物种净指数地图 | SI |
| K3 | fig_v4_species_trend_classification | 06 | 扩张/稳定/收缩柱状 | SI |
| K4 | fig_v4_species_trend_density | 06 | Theil-Sen 斜率密度 | SI |
| K5 | fig_v4_cwm_spatial_traits | 06 | CWM diet_specialization + habitat_breadth 地图 | SI |
| K6 | fig_v4_grid_size_sensitivity | 16 | 50 km vs 100 km 趋势配对 | SI |
| K7 | fig_v4_eps_sensitivity | 17 | 5 档 eps 相关性热力图 | SI |

---

## S6 Tables of all CSV outputs

按 `results_v4/` 表名 + 简短描述列表（约 30 张表），含：
- table_diversity_summary_<run_label>.csv
- table_trend_summary_<run_label>.csv
- table_baselga_summary_<run_label>.csv & table_baselga_global_*
- table_mann_kendall_*
- table_naive_vs_corrected_*
- table_species_trend_classify_*
- table_species_hotspot_*
- table_cwm_spatial_pattern_*
- table_brms_driver_coefs_*
- table_brms_driver_R2_*
- table_loo_driver_*
- table_trait_regression_coefs_*
- table_trait_regression_R2_*
- table_trait_regression_loo_comparison_*
- table_grid_size_sensitivity_50km_vs_100km
- table_eps_sensitivity_correlation
- table_eps_sensitivity_trend_direction
- table_breeding_filter_audit_v4
- table_convergence_diagnostics_*
- table_convergence_summary_*
- table_ppc_bayesian_pvalue_*
- table_ppc_global_summary_*

---

## S7 Software environment

- R version + sessionInfo()
- spOccupancy, brms, cmdstanr, loo, DHARMa, qs, sf, terra, ape, picante 版本
- Zenodo DOI（投稿前）
