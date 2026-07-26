# 中国鸟类 2000–2024 动态占域研究 — v4 全面审查与优化报告

> 审查时间：2026-05-11  审查人：生态学领域资深审稿人 / 高级编辑视角
> 项目根：`/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis`
> 审查对象：v1 / v2 / v3 三套代码与产出
> 配套交付：`code_v4/`、`results_v4/`、`figures_v4/`、`data/external/README_v4.md`

---

## 0 摘要（Executive Summary）

本项目立意优秀、方法学创新性强（首次在全国尺度上做"占域校正 + 后验传播 + 群落性状/系统发育多维度时空动态"），但**当前完成度与可投状态之间还有显著差距**。核心结论：

1. **方法学层**：v3 的设计方向正确，但**关键的不确定性传播链条尚未完整跑通**——FULL stMsPGOcc（200 物种 × 4 链）未完成，下游 05/06 因此没有 v3 版本结果与图集。已有 v2 结果可用作过渡，但**文稿数字与最新模型脱节**。
2. **数据层**：观鸟平台 + eBird 合并 + source-aware 去重是亮点；但**繁殖季过滤的时间偏差审计、100 km 尺度敏感性、西部空间覆盖缺失**这三处审稿人必问的问题缺失证据材料；外部栅格依赖（WorldClim/landcover/texture）部分缺失，导致 10 km 网格无法启用。
3. **代码层**：v3 解决了 v1/v2 的几处硬伤（PD 阈值化思路、log10 双峰、丰富度 bug、OLS 端点敏感），但**新引入 3 处 bug**：`lo()` 漏依赖（LOO 永远 NULL）、单链 R-hat 假报 1.0（误报"已收敛"）、`combine_chains()` 对 psi.samples 4D 合并未 flatten chain 维度（下游 05 维度错位风险）。
4. **图表层**：Nature 风格函数齐全，但 `figures_v3/` 完全为空；`standard_occupancy_figures.md` 列的 31 张标准图集 v3 至少缺 12 张（MCMC 诊断、PPC、物种系数 heatmap、AR1 ρ、NMDS+envfit 等）。
5. **文稿层**：v2 草稿骨架完整但**数字版本混杂**（886 候选物种 / 200 建模物种 / 60 pilot 三处共存），detection 协变量描述与实际公式不符，PD 描述与代码不符，讨论段薄、缺与已发表中国鸟类研究的对比，因果语言过强；缺英文 abstract / cover letter / supplementary outline。

整体投稿成熟度评估（按 Global Change Biology / Nature Communications 的尺度）：**60/100**。
按 P0→P3 路线推进，**6–8 周可达到 85/100** 的"可投"状态。

---

## 1 数据层（Data layer）

### 1.1 优点
- 双源（中国观鸟记录平台 + eBird/GBIF）合并 → 1448 万 → 749 万去重事件，**数据基础在国内同类研究中规模最大、覆盖最长**。
- v3 把 v2 的"匿名用户 anon_lon_lat_date 合并不同观测者"bug 修正为 source-aware hash 去重（`01_merge_birdwatch_ebird.R` L127–133）。
- 候选物种筛选阈值（n_grid_year_detections ≥ 80, n_blocks_detected ≥ 2, n_grids_detected ≥ 10）合理。

### 1.2 问题与建议（按优先级）

#### D1【P1，方法学要害】繁殖季过滤的时间偏差未审计
- **现象**：`02_build_survey_history.R` L36–37 强制保留 4–8 月。但中国观鸟平台早期（2000–2010）以候鸟越冬观测为主，过滤后 P1/P2 周期事件可能损失 >50%；P4/P5 几乎不掉，**人为放大时间趋势**。
- **审计脚本**：`02_build_survey_history.R` L54–76 已写 `breeding_audit`，但 `results_v3/table_breeding_filter_audit.csv` 不存在（02 未在 v3 环境下完整跑过）。
- **必做**：
  1. 跑 `code_v4/02_build_survey_history.R`（v4 版本继承 v3 但强制写出审计表），保存 `results_v4/table_breeding_filter_audit_v4.csv`；
  2. 执行 `code_v4/15b_sensitivity_breeding_season.R`：跑一版"全年 + month as detection covariate"对比模型，输出 `table_breeding_vs_year_round_trend_consistency.csv`；
  3. Methods 写入 "Breeding-season filter justification + audit"；Supplementary 放审计表。

#### D2【P2】西部 / 中部空间覆盖严重不足
- **现象**：1261 / 1739 网格有数据，缺失集中西藏、新疆、青海、四川山区。
- **建议**：
  - 所有趋势图把"无数据网格"显式标灰色（`na.value="grey92"`），不要让 `inner_join` 偷偷踢掉；
  - Methods + Discussion 写明 "inference scope: 1261 sampled grids; no extrapolation to unsampled western mountains"；
  - Supplementary 出 `fig_v4_sampling_coverage_audit_map.png`：5 主期 × 网格的 visit 数。

#### D3【P2】100 km 网格尺度过粗
- **现象**：00_config.R 已开放 `V3_GRID_SIZE_KM` 但 10 km / 50 km 从未运行；100 km 一格能横跨海拔 500–3000 m。
- **建议**：跑 `code_v4/16_sensitivity_grid_size.R`：top60 sp 在 50 km 网格做 pilot，输出 `table_grid_size_sensitivity_50km_vs_100km.csv`；Discussion 加 "coarse-grained landscape occupancy"。

#### D4【P2】检测协变量兼容性 has_duration 全 0 fallback
- **现象**：`03_prepare_environment.R` L141 当 v2 visit_effort overlap <30% 时直接全 0 fallback，导致 duration 信号丧失。
- **建议**：v4 中 `03_prepare_environment.R` 加入显式 warning + 输出 `table_detection_covariate_audit.csv`（按 period × grid 展示 duration 缺失率）。

#### D5【P1，可立刻修复】外部栅格依赖断裂
- **现象**：
  - `data/external/worldclim/` 不存在 → bio4/7/11/13 在 10 km 场景下回退到 NA；
  - `landcover_built / cropland`、`texture_shannon`、`habitat_diversity_shannon`、`elev_sd` 在新提取流程里全为 NA（03 L172–176）；
  - 邻近项目硬编码：`TRAIT_IMPUTED_PATH` 指向 `bird_grid_community_analysis/`（00_config.R L128）。
- **建议**：
  - 按 `data/external/README_v4.md` 清单补齐 WorldClim BIO、ESA WorldCover、EarthEnv texture、SRTM；
  - 把 trait 表本地化到 `data/external/traits/`（即把 `table_species_traits_imputed.csv` 复制过来），`code_v4/00_config.R` 改路径。

#### D6【P3】去重 v2 vs v3 差 647 条无审计说明
- 建议跑 `code_v4/01b_dedup_diff_audit.R` 输出 `table_dedup_v3_vs_v2_audit.csv` 说明差异来源。

---

## 2 代码层（Code layer）

### 2.1 优点
- 模块化清晰：utils_core / utils_paths / utils_diversity / utils_mapping / utils_spatial / utils_diagnostics / utils_importance / utils_plots_advanced 八个工具模块；
- v3 修复 v2 的 log10 双峰 bug（`utils_diversity.R::prepare_trait_matrix` L260–305）；
- v3 修复丰富度计算（`div_taxonomic` 改 `sum(psi)`）；
- v3 改用 cmdstanr 后端（`05_postprocess_diversity.R` L502）；
- v3 加入 Theil-Sen + Mann-Kendall + Baselga + 后验传播框架。

### 2.2 问题清单（含 v4 修复方案，文件落到 `code_v4/`）

| # | 严重度 | 文件:行 | 问题 | v4 修复方式 |
|---|---|---|---|---|
| C1 | **P1** | `code_v3/05_postprocess_diversity.R:522` & `14_species_trait_regression.R:157` | `lo(brms_fit)` 应为 `loo::loo(brms_fit)`；被 tryCatch 吞掉，**所有 LOO 永远是 NULL** | `code_v4/05_postprocess_diversity.R` 顶部加 `library(loo)`，调用 `loo::loo(brms_fit)`；写出 `table_loo_*.csv` |
| C2 | **P1** | `code_v3/04b_recover_diagnostics.R:432-460` | `compute_rhat()` 在单链时 `length(chain_list) < 2` 直接返回 1，**误报"已收敛"** | `code_v4/utils_diagnostics.R` 改为单链时返回 `NA` + warning |
| C3 | **P1** | `code_v3/04b_recover_diagnostics.R:207-214` | `abind_chains` 对 4D `psi.samples` 合并到 5 维但下游 05 用 `psi_samples[draw, sp, site, period]`，**维度错位风险** | `code_v4/04b_recover_diagnostics.R` 合并后立刻把 chain 维 flatten 进 draw 维（reshape `[draw, chain, ...]` 成 `[draw*chain, ...]`） |
| C4 | **P2** | `code_v3/utils_diversity.R:104-145` | `div_phylogenetic()` 实际用 `picante::pd(samp=1, ...)` × `sum(psi)`，**不是 manuscript 描述的概率加权 PD** | `code_v4/utils_diversity.R::div_phylogenetic_prob()` 按枝长 × 概率展开：`PD_prob = Σ_e L_e × (1 − ∏_{sp∈desc(e)}(1 − ψ_sp))` |
| C5 | **P2** | `code_v3/03_prepare_environment.R:83 vs 162` | `extract_env_from_raster()` 调用在 L83，定义在 L162（lucky bug：当前路径不一定触发） | `code_v4/utils_spatial.R::extract_env_from_raster()` 抽离到 utils |
| C6 | **P2** | `code_v3/05_postprocess_diversity.R:478-481` | brms 公式同时塞 11 个标准化协变量 + `gp(centroid_lon, centroid_lat, k=10)`，空间混淆 + 无变量选择 | `code_v4/05_postprocess_diversity.R` 改用 `horseshoe(df=1, par_ratio=0.1)` 稀疏先验；报告 marginal vs conditional R²；跑 k=10/20/50 三档 |
| C7 | **P1** | 04 输出 21.9 GB（pilot），FULL 估计 >80 GB | `psi.samples` + `z.samples` + `w.samples` 都保留 | `code_v4/04_run_stMsPGOcc_main.R`：(1) 在保存前显式 `fit$z.samples <- NULL; fit$w.samples <- NULL`；(2) 用 `qs::qsave(..., preset="fast")`；(3) 立刻 thin 一次保存 `psi_samples_thinned_v4_*.qs`，下游不再读完整 fit |
| C8 | **P3** | `code_v3/utils_diversity.R:524-579` | `theil_sen_slope` / `mann_kendall_test` O(n²) 纯 R 双循环 | `code_v4/utils_diversity.R` 改用 `mblm::mblm()` + `Kendall::MannKendall()`；保留纯 R 版本做 fallback |
| C9 | **P3** | 多处 utils 重复定义 `find_project_root` fallback | source-order 风险 | `code_v4/utils_paths.R::find_project_root()` 单点定义，其他 utils source 之 |
| C10 | **P3** | 候选物种 fallback 重复 | `04_run_stMsPGOcc_main.R:63-72`、`15_sensitivity_3yr_window.R:67-74` | 抽 `code_v4/utils_core.R::limit_candidate_species(max_n=200)` |
| C11 | **P3** | 种子不一致 | prepare 20260501、main 20260502；missForest/ranger/brms 无局部 set.seed | `code_v4/utils_seeds.R::set_seeds(stage)` 按 stage 派生 |
| C12 | **P2** | `eps = 1e-12` 在多样性函数里全局滥用 | 几乎所有物种被算"存在"，但 picante::pd 在二值层面失去占域校正含义 | C4 一并修；额外加 `code_v4/17_sensitivity_eps_threshold.R` 做 1e-12 / 1e-6 / 0.01 / 0.05 / 0.10 五档敏感性 |

### 2.3 v4 代码结构（与 v3 完全并行，不覆盖）

```
code_v4/
├─ 00_config.R                     # 路径本地化、driver 标签、Nature 规范
├─ utils_paths.R                   # find_project_root, v4_file, log_path
├─ utils_core.R                    # safe_read, qs_save_safe, limit_candidate_species, log_time
├─ utils_seeds.R                   # set_seeds(stage)
├─ utils_spatial.R                 # project_china_albers, extract_env_from_raster
├─ utils_diversity.R               # 修 PD 为真正概率加权；mblm/Kendall 替换
├─ utils_diagnostics.R             # 单链 R-hat 返回 NA；DHARMa gate
├─ utils_mapping.R                 # 沿用 v3，封 coord_sf bbox + 无十段线
├─ utils_plots_advanced.R          # 沿用 v3
├─ utils_importance.R              # 沿用 v3
├─ 01b_dedup_diff_audit.R          # 新增：v2 vs v3 差异审计
├─ 02_build_survey_history.R       # v3 基础 + 强制写出 breeding audit
├─ 03_prepare_environment.R        # 修 extract_env_from_raster 调用 + duration 审计
├─ 04_run_stMsPGOcc_main.R         # qs 保存 + thin + 内存优化
├─ 04b_recover_diagnostics.R       # psi.samples 4D chain flatten + 单链 R-hat fix
├─ 05_postprocess_diversity.R      # loo + horseshoe + CWM 空间表
├─ 05b_mcmc_diagnostic_plots.R     # 新增：R-hat/ESS/trace/density 4 类图
├─ 05c_ppc_bayesian_pvalue.R       # 新增：ppcOcc + fit.y vs fit.y.rep
├─ 06_figures_publication.R        # 加 species coef heatmap + AR1 ρ + CWM 空间
├─ 06b_regenerate_driver_plots.R   # 沿用 v3
├─ 07_render_manuscript.R          # 统一 v4 数字
├─ 08_render_journal_docx.R        # journal 提交版
├─ 14_species_trait_regression.R   # loo + LOO 比较有/无 phylo
├─ 15_sensitivity_3yr_window.R     # 沿用 v3
├─ 15b_sensitivity_breeding_season.R  # 沿用 v3 但加 strict 对比
├─ 16_sensitivity_grid_size.R      # 新增：50 km 网格 pilot 对比
├─ 17_sensitivity_eps_threshold.R  # 新增：eps 5 档敏感性
└─ README_v4.md                     # 执行顺序 + 验证清单
```

---

## 3 图表层（Figure layer）

### 3.1 优点
- `code_v3/utils_mapping.R` 已封装 Nature 风格（Arial 7.5 pt、89/120/183 mm、scale_fill_nature_c、theme_nature_map）；
- v2 已有部分图通过 `map_quality_rules.md` 7 条硬规则（无十段线、coord_sf 固定 bbox、无 NA 子图、color scale within-metric z-score）；
- 雨林图 / 蜂群图工具已在 `utils_plots_advanced.R`。

### 3.2 问题与缺失

| 图组 | 标准图（来自 `standard_occupancy_figures.md`） | v3 状态 | v4 补齐路径 |
|---|---|---|---|
| **A. MCMC 收敛** | R-hat 直方、ESS 直方、trace、posterior density 4 链叠合 | 缺 | `code_v4/05b_mcmc_diagnostic_plots.R` 全部出 |
| **B. PPC** | Bayesian p-value（Freeman-Tukey + Chi²）、fit.y vs fit.y.rep 散点 | 缺 | `code_v4/05c_ppc_bayesian_pvalue.R` |
| **C. 群落层** | beta.comm / alpha.comm caterpillar、hierarchical shrinkage | 部分（CSV 有，图无） | `code_v4/06_figures_publication.R` 加 caterpillar |
| **D. 物种层** | 物种 × 协变量 heatmap，× 标 95% CRI 不跨 0 | **缺** | `code_v4/06_figures_publication.R` 新加 H 段 |
| **E. 时空动态** | 物种 psi(t) 轨迹、群落 psi/det 时间切片地图、naive vs corrected 散点 | v3 06 已写 D 段，但因 FULL 缺失未跑 | FULL 跑完即出 |
| **F. AR1 时序** | rho 后验密度（按物种）、sigma.sq.t 后验 | 缺图 | `code_v4/06_figures_publication.R` 加 I 段 |
| **G. 多样性传播** | richness/Shannon/PD/Rao Q 时序带 95% CRI 阴影、趋势地图 + CRI 宽度 | 06 E 段已写 | FULL 跑完即出 |
| **H. 群落重组** | Bray temporal beta raincloud、Baselga turnover vs nestedness 双面板地图、NMDS + envfit FDR | 缺 NMDS | `code_v4/06_figures_publication.R` 加 J 段 |
| **I. 驱动回归** | brms 系数 forest plot、conditional_effects、DHARMa 4 panel、残差空间自相关地图 | 06b 已写、DHARMa 已写 | 验证 4 panel 完整 |
| **J. 性状驱动** | 性状 ~ trend forest plot、PGLS / brms phylogenetic 残差 | 缺图 | `code_v4/14_species_trait_regression.R` 后段加出图 |

### 3.3 地图自检清单（7 条硬规则，每张图都过一遍）
1. `coord_sf(xlim=c(73,135), ylim=c(18,54), expand=FALSE)` 固定 bbox
2. 不画十段线 / 鹰眼图 / inset
3. `facet_wrap` 无 NA 子图（`metric` 转 factor + 指定 levels）
4. 缺值网格保留 NA → `na.value="grey92"`，不被 inner_join 踢掉
5. 跨指标 facet 用 within-metric z-score（clip ±2.5）；非负量用 `scico::lajolla`
6. 南方（Yunnan/Guangxi/Hainan/Taiwan/Guangdong/Fujian）有数据网格不要空白
7. 北方（Heilongjiang/Inner Mongolia/Xinjiang）无横向条纹

---

## 4 文稿层（Manuscript layer，编辑视角）

### 4.1 v2 草稿（`results_v2/manuscript_v2_zh.md`）评分

| 维度 | 评分 | 问题 |
|---|---|---|
| 引言 narrative | B− | 科学问题 Q1–Q5 清晰，但缺"已有方法的局限 → 本研究的突破点"对比叙事 |
| Methods 完整度 | C+ | 数字 886 / 200 / 60 三个版本混杂；detection 公式与 v3 实际不符；PD 描述与代码不符；MCMC 配置与 pilot 不符 |
| Results 数字 | B− | 5 个 period 的 richness / PD / Shannon 给得齐，但缺驱动因子结果、缺物种层（Q4） |
| Discussion | D | 仅 3 句话，无文献对比，因果语言过强 |
| Limitations | C | 列了 3 条但浅；缺 closure assumption、采样偏差、外部验证 3 条 |
| 图表呼应 | C | 正文未引用任何 figure / table 编号 |
| 英文摘要 / cover letter / SI | F | 全无 |

### 4.2 v4 修订要点（落到 `results_v4/manuscript_v4_zh.md` + 4 个新文件）

#### M1【P0】统一数字版本
- 删除所有 "886 候选物种" 出现的结果段；
- Methods 明确两步筛选：(1) 候选池 886 通过 detection threshold；(2) 建模 200 取 top reach-rate；
- 给出 200 vs 全 886 的敏感性（在 `15_sensitivity_3yr_window.R` 之外新增 candidate threshold sensitivity）。

#### M2【P1】PD 描述与代码一致
- 选项 A（推荐）：v4 把 `div_phylogenetic_prob()` 改为真正概率加权（C4 修复）；
- Methods 公式与算法描述完全对齐。

#### M3【P1】Detection 协变量描述更新
- Methods：`p ~ log_events + log_duration + has_duration`（`has_duration` 作为 missingness indicator）。

#### M4【P1】MCMC 配置描述更新（以 FULL 实际为准）
- `4 chains × 400 batch × 25 length, burn=5000, thin=2`；
- 收敛标准 R-hat ≤ 1.05, ESS ≥ 200；FULL 跑完后写实际百分比。

#### M5【P2】Discussion 扩展至 5 段
1. 占域校正 vs 朴素方法的量化差异（对比 BBS 类研究：Yang et al. 2020 Diversity & Distributions; Liang et al. 2018 Biol Conserv）；
2. Turnover vs Nestedness 比例及与全球 meta-analysis 对比（Baselga & Orme 2012 MEE; Soininen et al. 2018 Ecography）；
3. 东西部气候 vs 人为驱动差异；HFI 变化与建成度对边缘物种的影响（参 Tucker et al. 2018 Science）；
4. 性状（HWI、diet_specialization、habitat_breadth）的扩散/特化解释力 vs Pigot et al. 2020 Nat Ecol Evol；
5. 局限：closure assumption、采样偏差、缺乏外部验证、相关性 ≠ 因果。

#### M6【P2】因果语言修订
- 全文 `drive / cause` → `associate with / explain spatial variation in`；
- 加 "We acknowledge that our analysis identifies correlations, not causal effects, due to the observational nature of citizen-science data."

#### M7【P3】标题更新
建议候选（任选一）：
- "Citizen-science reveals occupancy-corrected biotic homogenization of Chinese birds in a warming century (2000–2024)"
- "Detection-corrected community dynamics expose hidden range shifts and functional reshuffling of Chinese avifauna"

#### M8【P1】英文 abstract + cover letter + SI outline
- `results_v4/manuscript_v4_en_abstract.md`：250 词 structured；
- `results_v4/cover_letter_v4.md`：1 页，强调"首次全国尺度多十年占域校正群落动态"；
- `results_v4/supplementary_outline_v4.md`：S1 数据流，S2 MCMC 诊断（A 组），S3 PPC（B 组），S4 敏感性（3yr / breeding / 50km / eps），S5 完整图集（31 张），S6 表格清单。

---

## 5 优先级路线图（按 4–8 周节奏）

### Week 1（P0：把 pipeline 跑通）
- [ ] 用 `code_v4/04_run_stMsPGOcc_main.R` 在 1 TB 服务器（<SERVER_IP>）上跑 FULL 200 sp × 4 chain；
- [ ] `code_v4/04b_recover_diagnostics.R` 合并 + 收敛；要求 max R-hat ≤ 1.05；
- [ ] `code_v4/05_postprocess_diversity.R` 出 diversity / Baselga / trend / species_trend / naive_vs_corrected / hotspot；
- [ ] `code_v4/06_figures_publication.R` 出第一批 v4 图集；
- [ ] `code_v4/07_render_manuscript.R` 渲染 v4 中文稿；
- [ ] 修 lo()→loo::loo()、单链 R-hat→NA、psi.samples chain flatten 三个 bug。

### Week 2（P1：诊断 + 敏感性）
- [ ] `code_v4/05b` MCMC 4 类诊断图；`code_v4/05c` PPC；
- [ ] `code_v4/15b` 繁殖季敏感性；`code_v4/16` 50 km 网格 pilot；`code_v4/17` eps 敏感性；
- [ ] `code_v4/14` 性状回归含/不含 phylo random effect 的 LOO 比较；
- [ ] 物种 × 协变量 heatmap、AR1 ρ 后验密度图。

### Week 3（P1：稿件迭代）
- [ ] 数字版本统一（M1）；PD 描述对齐（M2）；detection 公式更新（M3）；MCMC 配置更新（M4）；
- [ ] 英文摘要 + cover letter + SI outline；
- [ ] 图集自检（7 条硬规则）。

### Week 4–6（P2：鲁棒性 + 论文打磨）
- [ ] Discussion 5 段扩展；
- [ ] NMDS + envfit (FDR)；
- [ ] brms horseshoe 稀疏先验 + GP k=10/20/50 比较；
- [ ] varpart vs RF 互补图；
- [ ] 31 张标准图 checklist 全勾完；
- [ ] 与已发表中国鸟类研究 3–5 区域的趋势比较表。

### Week 7–8（投稿准备）
- [ ] 外审风格 code review（请 1 位 co-author 跑一遍 `code_v4/01..17.R`）；
- [ ] cover letter 抛光；
- [ ] 期刊格式（Word + figures + tables）；
- [ ] GitHub 仓库 + Zenodo DOI + 数据 README。

---

## 6 风险与降级方案

- **风险 A**：FULL stMsPGOcc 在 24 GB 本机 OOM。
  **降级**：先 `top100sp_ar1`，再扩 200；或先 `tMsPGOcc`（无空间）跑通整套 pipeline 再加空间项。
- **风险 B**：brms + GP 在 1308 网格上编译/采样 >24 h。
  **降级**：换 `INLA` 或 `mgcv::gam(... + s(lon, lat, bs="gp"))` 做快诊；brms 只用于核心 4 个指标。
- **风险 C**：审稿人质疑 closure assumption。
  **应对**：D1 繁殖季 + 3 yr 窗口 + month-as-covariate 三个 robustness 摆 SI，主文一句话引用。
- **风险 D**：FULL 跑完 R-hat 仍 >1.05。
  **应对**：增加 batch 到 600、burn 到 8000；或降低 n.factors=3 减少 latent 维度。

---

## 7 验证（v4 完成的判定标准）

1. **MCMC**：`table_convergence_diagnostics_v4_full_*.csv` 中 max R-hat ≤ 1.05，min ESS ≥ 200。
2. **PPC**：Freeman-Tukey Bayesian p-value 0.1–0.9。
3. **后处理表齐全**（`results_v4/`）：
   - `table_diversity_summary_v4_full_*.csv`
   - `table_baselga_summary_v4_full_*.csv` + `table_baselga_global_*`
   - `table_trend_summary_v4_full_*.csv`（含 ols + theil_sen 两列）
   - `table_mann_kendall_v4_*.csv`
   - `table_species_trend_classify_v4_*.csv`
   - `table_naive_vs_corrected_v4_*.csv`
   - `table_loo_*.csv`（确认非 NULL）
4. **图集齐全**（`figures_v4/`）：≥ 15 张，覆盖 A–J 10 类至少 8 类。
5. **地图 7 条规则全过**。
6. **稿件数字一致性**：`grep -oE "[0-9,.]+%?" results_v4/manuscript_v4_zh.md` 全部能在 `table_*` CSV 对齐。
7. **敏感性结果**：3 yr / 繁殖季 / 50 km / eps 四个敏感性的趋势方向一致性 ≥ 70%。
8. **可复现性**：干净环境下 `Rscript code_v4/01..17.R` 全跑通无报错。

---

## 8 交付清单（本次 v4 落地）

本审查附带产出以下文件（全部在 v4 独立路径，不动 v3）：

```
results_v4/
├─ REVIEW_v4_full_audit_cn.md       ← 本文件（综合审查报告）
├─ manuscript_v4_zh.md              ← 中文稿件骨架（统一数字 + 公式 + 文献）
├─ manuscript_v4_en_abstract.md     ← 英文 250-word structured abstract
├─ cover_letter_v4.md               ← 1 页 cover letter
└─ supplementary_outline_v4.md      ← SI S1–S6 大纲

code_v4/
├─ README_v4.md                     ← 执行顺序 + 验证清单
├─ 00_config.R                       ← 本地化路径 + Nature 规范
├─ utils_paths.R / utils_core.R / utils_seeds.R / utils_spatial.R
├─ utils_diversity.R                 ← 概率加权 PD + mblm/Kendall
├─ utils_diagnostics.R               ← 单链 R-hat → NA
├─ 04_run_stMsPGOcc_main.R           ← qs + thin + 内存优化
├─ 04b_recover_diagnostics.R         ← psi.samples chain flatten
├─ 05_postprocess_diversity.R        ← loo::loo + horseshoe + CWM 空间
├─ 05b_mcmc_diagnostic_plots.R       ← MCMC 4 类诊断图
├─ 05c_ppc_bayesian_pvalue.R         ← ppcOcc + fit.y vs rep
├─ 06_figures_publication.R          ← 含 species coef heatmap + AR1 ρ + CWM
├─ 14_species_trait_regression.R     ← loo 比较 phylo random
├─ 16_sensitivity_grid_size.R        ← 50 km pilot
├─ 17_sensitivity_eps_threshold.R    ← 5 档 eps 敏感性
└─ 07_render_manuscript.R + 08_render_journal_docx.R

data/
└─ external/
   └─ README_v4.md                   ← 外部栅格/性状/系统发育数据来源 + DOI

figures_v4/
└─ （等 FULL 跑完后由 code_v4/06 出图填充）
```

---

> **下一步**：建议先看 `code_v4/README_v4.md`（执行顺序）+ `results_v4/manuscript_v4_zh.md`（新数字版本）；然后把 `code_v4/04_run_stMsPGOcc_main.R` 同步到服务器跑 FULL；FULL 落盘后回到本机依次跑 04b → 05 → 05b → 05c → 06 → 07，约 1 天即可拿到 v4 全套图集 + 稿件骨架。
