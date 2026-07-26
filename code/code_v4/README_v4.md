# code_v4 — 执行顺序与验证清单

> v4 是 v3 的"全面审稿对齐 + bug 修复"版本，与 v3 完全并行（不覆盖任何 v3 文件）。
> 所有 v4 输出固定到 `code_v4 / results_v4 / figures_v4 / data/derived_v4 / logs_v4`。
> 配套交付物 `results_v4/REVIEW_v4_full_audit_cn.md` 解释每条修复对应的 v3 问题编号。

---

## 0. 环境准备

### R 包
```r
# 主流水线
install.packages(c("spOccupancy", "brms", "cmdstanr", "loo",
                    "DHARMa", "abind", "qs",
                    "sf", "terra", "geosphere",
                    "dplyr", "tidyr", "tibble", "readr", "stringr",
                    "ggplot2", "patchwork", "scales",
                    "ape", "picante", "vegan", "Kendall", "mblm",
                    "scico", "digest", "glue", "forcats"))
# cmdstanr 的 cmdstan 实际安装
cmdstanr::install_cmdstan()
```

### 外部数据
按 `data/external/README_v4.md` 检查清单补齐外部栅格与性状文件。

### 环境变量
```bash
export BIRD_PROJECT_ROOT="/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
export V4_CODE_DIR="$BIRD_PROJECT_ROOT/code_v4"

# 可选：覆盖默认参数
export V4_GRID_SIZE_KM=100      # 主跑 100 km；敏感性时改 50
export V4_MAX_SPECIES=200       # full run 物种数
export V4_PILOT=0                # 1 = pilot 20 sp（本机可跑）
export V4_N_OMP_THREADS=1
export V4_CHAIN_ID=              # 留空 = 单进程跑所有链；1..4 = 链并行
```

---

## 1. 执行顺序（Pilot 验证 → Full 跑通）

### 1.1 Pilot（本机 24 GB 即可，1 小时内）
```bash
cd "$BIRD_PROJECT_ROOT"
V4_PILOT=1 Rscript code_v4/01_merge_birdwatch_ebird.R  # 如有 v3 已合并文件可跳过
V4_PILOT=1 Rscript code_v4/02_build_survey_history.R
V4_PILOT=1 Rscript code_v4/03_prepare_environment.R
V4_PILOT=1 Rscript code_v4/04_run_stMsPGOcc_main.R     # pilot 20 sp × 4 链
V4_PILOT=1 Rscript code_v4/04b_recover_diagnostics.R
V4_PILOT=1 Rscript code_v4/05_postprocess_diversity.R
V4_PILOT=1 Rscript code_v4/05b_mcmc_diagnostic_plots.R
V4_PILOT=1 Rscript code_v4/05c_ppc_bayesian_pvalue.R
V4_PILOT=1 Rscript code_v4/06_figures_publication.R
V4_PILOT=1 Rscript code_v4/14_species_trait_regression.R
V4_PILOT=1 Rscript code_v4/07_render_manuscript.R
V4_PILOT=1 Rscript code_v4/08_render_journal_docx.R
```

### 1.2 Full（建议服务器 ≥ 256 GB；链并行）
```bash
# 在服务器（<SERVER_IP>）上：
# 链并行：4 个独立进程，每个 V4_CHAIN_ID=1..4
for i in 1 2 3 4; do
  V4_CHAIN_ID=$i Rscript code_v4/04_run_stMsPGOcc_main.R > logs_v4/04_chain${i}.log 2>&1 &
done
wait

# 合并 + 诊断
Rscript code_v4/04b_recover_diagnostics.R

# 下游
Rscript code_v4/05_postprocess_diversity.R
Rscript code_v4/05b_mcmc_diagnostic_plots.R
Rscript code_v4/05c_ppc_bayesian_pvalue.R
Rscript code_v4/06_figures_publication.R
Rscript code_v4/14_species_trait_regression.R
Rscript code_v4/07_render_manuscript.R
Rscript code_v4/08_render_journal_docx.R
```

### 1.3 敏感性分析（可选，但投稿前必跑）
```bash
# 3 年窗口
Rscript code_v4/15_sensitivity_3yr_window.R

# 繁殖季 vs 全年
Rscript code_v4/15b_sensitivity_breeding_season.R

# 50 km 网格（需先以 V4_GRID_SIZE_KM=50 跑 pilot 04）
V4_GRID_SIZE_KM=50 V4_PILOT=1 Rscript code_v4/02_build_survey_history.R
V4_GRID_SIZE_KM=50 V4_PILOT=1 Rscript code_v4/03_prepare_environment.R
V4_GRID_SIZE_KM=50 V4_PILOT=1 Rscript code_v4/04_run_stMsPGOcc_main.R
V4_GRID_SIZE_KM=50 V4_PILOT=1 Rscript code_v4/04b_recover_diagnostics.R
V4_GRID_SIZE_KM=50 V4_PILOT=1 Rscript code_v4/05_postprocess_diversity.R
Rscript code_v4/16_sensitivity_grid_size.R

# eps 阈值敏感性
Rscript code_v4/17_sensitivity_eps_threshold.R
```

---

## 2. v4 相对 v3 的关键修复

| ID | v3 问题 | v4 修复 | 文件 |
|---|---|---|---|
| C1 | `lo(brms_fit)` 漏依赖，所有 LOO=NULL | 顶部 `library(loo)`；调用 `loo::loo()` | 05, 14 |
| C2 | 单链 R-hat 假报 1.0 | 单链返回 NA + warning | utils_diagnostics |
| C3 | psi.samples 4D chain 未 flatten | 合并后立即 chain → draw flatten | 04b（flatten_chain_dim） |
| C4 | div_phylogenetic 用 picante::pd 二值化 | 真正概率加权 PD：PD = Σ_e L_e × (1 − ∏(1 − ψ)) | utils_diversity::div_phylogenetic_prob |
| C5 | extract_env_from_raster 定义在调用之后 | 抽到 utils_spatial.R | utils_spatial |
| C6 | brms 11 协变量 + GP 空间混淆 | horseshoe 稀疏先验 + GP k 三档（10/20/50）+ 报告 marginal vs conditional R² | 05 |
| C7 | psi.samples 22 GB 内存压力 | 保存前 z/w slots = NULL；qs::qsave preset="fast"；立刻 thin 落盘 | 04 |
| C8 | Theil-Sen / MK 纯 R O(n²) | mblm / Kendall 优先，纯 R fallback | utils_diversity |
| C9 | 性状回归不比较有/无 phylo | 同时拟合 M0/M1 + loo_compare | 14 |
| C10 | 候选物种 fallback 重复实现 | utils_core::limit_candidate_species 单点定义 | utils_core, 04, 15 |
| C11 | 种子不一致 | utils_seeds::set_seeds(stage) 按 stage 派生 | utils_seeds, 所有 stage |
| C12 | eps=1e-12 全局滥用 | PSI_EPS_DEFAULT=0.05 默认 + 17 做 5 档敏感性 | 00_config, utils_diversity, 17 |
| D5 | 邻近项目硬编码 | 00_config.R trait 路径本地化优先；data/external/README_v4.md | 00_config |
| M1-M4 | 稿件数字 / detection / PD 描述 / MCMC 不一致 | 07_render_manuscript.R 全部从 CSV 注入 | 07, manuscript_v4_zh.md |

---

## 3. 验证清单（投稿前逐项过）

### 3.1 MCMC
```bash
ls results_v4/table_convergence_diagnostics_v4_full_*.csv
# 检查：max R-hat ≤ 1.05；min ESS ≥ 200
```

### 3.2 PPC
```bash
ls results_v4/table_ppc_global_summary_v4_full_*.csv
# 检查：Freeman-Tukey median bp 在 0.1–0.9 之间
```

### 3.3 后处理产物完整
```bash
ls results_v4/ | grep -E "v4_full_200sp_ar1_spatial"
# 至少应有以下 11 张表：
# - table_diversity_summary
# - table_community_metrics_with_cri
# - table_temporal_dynamics_summary
# - table_baselga_summary & table_baselga_global
# - table_trend_summary
# - table_mann_kendall
# - table_naive_vs_corrected
# - table_species_trend & table_species_trend_classify
# - table_species_hotspot
# - table_cwm_spatial_pattern
# - table_brms_driver_coefs & table_brms_driver_R2 & table_loo_driver
# - table_trait_regression_coefs & table_trait_regression_R2 & table_trait_regression_loo_comparison
```

### 3.4 图集
```bash
ls figures_v4/ | grep "v4_" | wc -l
# 应该 ≥ 15，覆盖 A–K 11 组：
# A: mcmc_rhat_histogram, mcmc_ess_histogram, mcmc_trace_core_params, mcmc_density_core_params
# B: ppc_fit_y_vs_rep
# D: species_covariate_heatmap
# E: multidiversity_timeslices, multidiversity_uncertainty, multidiversity_trajectories
# F: mcmc_ar1_rho_density_by_species, cwm_spatial_traits
# G: community_trends_with_cri
# H: temporal_beta_baselga, baselga_proportion
# I: brms_driver_forest
# J: dharma_driver_* / dharma_trait_regression_*
# K: naive_vs_corrected_trend, species_hotspot_map, species_trend_classification, species_trend_density
#    grid_size_sensitivity, eps_sensitivity
```

### 3.5 地图 7 条硬规则（手工核对）
对每张地图：
1. 固定 bbox：xlim=c(73,135), ylim=c(18,54), expand=FALSE
2. 无十段线 / 鹰眼图 / inset
3. facet_wrap 无 NA 子图（metric 转 factor + 指定 levels）
4. 缺值 grid 显式灰色（na.value="grey92"），不被 inner_join 踢掉
5. 跨指标 facet 用 within-metric z-score（clip ±2.5）；非负量用 scico::lajolla
6. 南方（云南/广西/海南/台湾/广东/福建）有数据网格不空白
7. 北方（黑龙江/内蒙/新疆）无横向条纹

### 3.6 稿件数字一致性
```bash
grep -oE "[0-9,.]+%?" results_v4/manuscript_v4_zh.md | sort -u
# 每个数字应能在 results_v4/table_*.csv 中找到对应来源
```

### 3.7 敏感性结果一致性
```bash
ls results_v4/table_sensitivity_* results_v4/table_grid_size_sensitivity_* results_v4/table_eps_sensitivity_*
# 趋势方向一致性 ≥ 70%
```

### 3.8 可复现性
```bash
# 在干净 R session：
Rscript -e 'source("code_v4/00_config.R"); source("code_v4/utils_paths.R"); source("code_v4/utils_core.R"); source("code_v4/utils_diversity.R"); source("code_v4/utils_diagnostics.R"); cat("OK\n")'
# 无报错 → utils 链路通畅
```

---

## 4. 风险与降级

- **OOM (out of memory) on 04**：先用 `V4_PILOT=1 V4_MAX_SPECIES=60` 跑通；再扩到 100、200。
- **stMsPGOcc segfault**：降级到 tMsPGOcc（无空间），先验证数据管线；空间项后补。
- **brms + GP 编译 >24 h**：降到 k=10；或换 `mgcv::gam(... + s(centroid_lon, centroid_lat, bs="gp"))` 做快速诊断。
- **FULL R-hat > 1.05**：增加 batch 到 600、burn 到 8000；或降 n.factors=3。

---

## 5. 输出落地清单

```
data/derived_v4/
├─ survey_history_v4.rds
├─ grid_environment_v4.rds
├─ detection_covariates_v4.rds
├─ trait_extended_v4.rds                          # 来自 03b（待写或复制 v3）
├─ phylogeny_matched_v4.rds                       # 来自 v3/v2 复制
├─ stMsPGOcc_fit_<run_label>_chain1..4.qs         # 链并行模式
├─ stMsPGOcc_fit_<run_label>_combined.qs          # 04b 合并
├─ psi_samples_thinned_<run_label>_chain1..4.qs   # 04 立刻 thin
├─ psi_samples_thinned_<run_label>_combined.qs    # 04b 合并 thin
├─ brms_driver_*_k{10,20,50}_<run_label>.qs
├─ brms_trait_regression_M{0,1}_<run_label>.qs
├─ loo_*_<run_label>.qs
└─ varpart_*_<run_label>.qs

results_v4/
├─ REVIEW_v4_full_audit_cn.md                    # 本审查报告
├─ manuscript_v4_zh.md / .docx
├─ manuscript_v4_en_abstract.md / .docx
├─ cover_letter_v4.md / .docx
├─ supplementary_outline_v4.md / .docx
├─ workflow_v4_zh.md
├─ table_model_summary_<run_label>.csv
├─ table_model_io_audit_*.csv
├─ table_beta_community / alpha_community / spatial_params / beta_species
├─ table_convergence_diagnostics / summary
├─ table_ppc_bayesian_pvalue / global_summary
├─ table_diversity_summary / community_metrics_with_cri
├─ table_temporal_dynamics_summary
├─ table_baselga_summary / baselga_global
├─ table_trend_summary / mann_kendall
├─ table_naive_vs_corrected
├─ table_species_trend / classify / hotspot
├─ table_cwm_spatial_pattern
├─ table_brms_driver_coefs / R2 / loo_driver
├─ table_trait_regression_coefs / R2 / loo / loo_comparison
├─ table_sensitivity_3yr_vs_5yr
├─ table_grid_size_sensitivity_50km_vs_100km
└─ table_eps_sensitivity_correlation / trend_direction

figures_v4/
└─ fig_v4_*_<run_label>.{png,pdf}                  # 共 15+ 张图，覆盖 A–K
```

---

## 6. 帮助与故障排查

如某 stage 报错：
1. 看 `logs_v4/<stage>.log` 末尾 50 行；
2. 检查 `results_v4/table_model_io_audit_*.csv` 确认输入路径；
3. 对照本 README §2 表查找是否是 v3 → v4 修复中遗漏的边角；
4. 如需协助，把 `logs_v4/` 全部打包 + 报错具体行号发给 Chenchen。
