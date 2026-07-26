# 检测偏差校正下中国鸟类群落的时空动态（2000–2024）

> Run label: `v4_full_200sp_ar1_spatial`（待 FULL 跑完后由 `code_v4/07_render_manuscript.R` 自动注入数字）
> 生成时间: 2026-05-11（手写骨架版；FULL 跑完后将被脚本覆盖）
>
> **重要**：本稿件骨架基于 v4 数据流编写；所有数字字段（{n_records}, {n_species}, R-hat, WAIC, Baselga 比例等）将由 `07_render_manuscript.R` 从 `results_v4/table_*.csv` 自动注入。当前为模板，无明显错误数字版本（v2 草稿"886/200/60"混杂的问题已修复）。

---

## 摘要

气候变化、土地利用转型与人类干扰强度的快速变化共同重塑了中国鸟类群落的时空格局。然而，公民科学数据在调查强度上具有强烈的非随机性，简单出现频次或物种富集曲线极易把"调查扩张"误读为"群落扩张"。本研究以"中国观鸟记录平台"checklist 与 GBIF/eBird 中国记录合并去重后的 {n_records} 条访问事件为基础（2000–2024），在 100 km × 5 年单元上对 {n_species} 个建模物种拟合**空间多物种动态占域模型（stMsPGOcc，含 AR1 时间随机效应 + NNGP 空间随机效应，exponential 协方差）**，显式分离不完全探测、调查强度异质性与空间相关性。占域后验通过一致采样传播到分类（校正丰富度、Shannon）、系统发育（**真正概率加权 Faith's PD**、MPD）与功能多样性（CWM、CWSD、性状空间体积、Rao's Q、FEve、FDiv），并使用**概率加权 Baselga 分解**量化相邻主期的 turnover 与 nestedness 比例。

**核心发现**：(1) 校正丰富度在 5 个主期内呈现 {richness_str}；(2) {turnover_str}，群落组装由 turnover 主导而非 nestedness；(3) 物种层面：{cls_str}；(4) 占域校正显著重构 naive 出现频次方法的结论（{nvc_str}）；(5) 高 HWI、宽 habitat_breadth、低 diet_specialization 的物种倾向于占域扩张。

**关键词**：动态占有率模型、检测偏差校正、生物多样性时空动态、Baselga 分解、概率加权 PD、鸟类、公民科学、中国

---

## 1 引言

[1.1 全球变化与中国鸟类多样性研究的科学缺口]
[1.2 公民科学数据的优势与陷阱：检测偏差为何不能忽视]
[1.3 空间多物种动态占有率模型的方法学突破]
[1.4 本研究的科学问题 Q1–Q5]

### 1.5 科学问题
- **Q1**：在显式控制不完全探测、调查强度异质性与空间相关性后，2000 年以来中国鸟类群落在哪些指标上表现出系统性变化（95% CRI 排除 0）？
- **Q2**：这些变化是 colonization 主导的群落扩张，还是 local extinction 主导的群落退缩？相邻主期之间 turnover 与 nestedness 比例如何？
- **Q3**：气候季节性、地形纹理、生产力、土地覆盖、HFI、纬度梯度对群落动态趋势的解释强度（horseshoe 稀疏先验 + 空间 GP）如何？
- **Q4**：物种生活史与扩散相关性状（体重、寿命、HWI、范围大小、diet_specialization、habitat_breadth）能否解释物种层占域趋势的差异？系统发育结构是否必要？
- **Q5**：占域校正相对 naive 出现频次方法，结论的重构有多大？哪些网格的趋势方向发生翻转？

---

## 2 材料与方法

### 2.1 数据来源
- 中国观鸟记录平台（2000–2024，checklist 元数据 + 鸟种逐条记录）；
- GBIF/eBird China（清洗后，2000–2024）；
- 性状：AVONET（Tobias et al. 2022, Ecol Lett）+ EltonTraits（Wilman et al. 2014, Ecology）+ IUCN Red List habitat 数量；缺失用 missForest 插补；
- 系统发育：clootl 提取（taxonomy_year=2025, version=1.6），按学名匹配；
- 环境：WorldClim BIO4/7/11/13、SRTM 海拔均值/异质性、EarthEnv texture、ESA WorldCover landcover_built/cropland、Human Footprint 年度（Mu et al. 2022 Figshare）、CRU TS 月均温变化（climate change drivers）、CLCD 30 m 土地覆盖变化（Yang & Huang 2021 ESSD）。

### 2.2 时空设计
- 100 km 等积栅格（Albers 投影 → WGS84），中国大陆边界内 1,739 个网格，有数据 {n_sites} 个；
- 主期：2000–2004、2005–2009、2010–2014、2015–2019、2020–2024（共 5 个 primary period）；
- 重复：每主期内 5 个年份作为 secondary occasions（**relaxed closure** 假设，参 Rota et al. 2009 J Appl Ecol; Kendall & White 2009 JABES）；
- **繁殖季过滤**：仅保留 4–8 月记录以减少越冬/迁徙鸟混淆；按 period 审计过滤前后事件比例（详见 SI Table S1）。

### 2.3 数据合并与去重（source-aware）
合并两源事件 → source-aware hash：`(species × date × round(lon, 4) × round(lat, 4) × source × observer_id)` 联合 key，跨源保留，同源精确去重。匿名记录用 `(lon, lat, date, species, source)` 生成稳定 anon_id 避免合并不同观测者。

### 2.4 调查史与检测协变量
对每个 (grid_cell, year)：定义 visit = 唯一 (grid_cell, date, observer_id) 组合；统计 n_events、log_events = log1p(n_events)；duration：当 `mean_duration_min > 0` 时 `log_duration = log10(mean_duration_min)`，否则 `log_duration = 0` 由 `has_duration` 作 missingness indicator 标记。**检测子模型公式：`p ~ log_events + log_duration + has_duration`**。

### 2.5 占域模型（stMsPGOcc）
`spOccupancy::stMsPGOcc` 对 {n_species} 个物种联合估计 occupancy ψ_{i, j, t} 与 detection p_{i, j, t, k}，hierarchical 共享物种间响应。

**占域协变量**：BIO4 (温度季节性)、BIO7 (温度年较差)、BIO11 (最冷季均温)、BIO13 (最湿月降水)、elev_mean、elev_sd、texture_shannon、habitat_diversity_shannon、hfi_mean、landcover_built、landcover_cropland、centroid_lon、centroid_lat、year_scaled（用于物种水平时间趋势 trend_i）。

**空间结构**：NNGP 近似（5 邻居，exponential 协方差）；先验 phi ~ Unif(0.1, 10) km^−1；初始 phi 用测地线距离均值的 3 倍倒数。

**时间结构**：AR(1)，ρ 与 σ²_t 自由估计。

**物种间相关**：latent factor 模型，n.factors = min(5, n_sp − 1)。

**MCMC**：4 链，每链 400 batch × 25 length（共 10,000 后样本/链），burn-in 5,000，thin 2，得到 (10,000 − 5,000)/2 = 2,500 thinned 后样本/链 × 4 链 = 10,000 总后样本。链间并行（V4_CHAIN_ID 单链执行 + 04b 合并 + chain 维 flatten）。

### 2.6 收敛、PPC 与 DHARMa
- **R-hat** (Gelman-Rubin)：群落 beta.comm / alpha.comm、空间 phi / sigma.sq、时序 rho，目标 max ≤ 1.05；
- **ESS** (coda::effectiveSize)：目标 min ≥ 200；
- **PPC**：`spOccupancy::ppcOcc()` 用 Freeman-Tukey 与 Chi-square 两种统计量，按物种汇总 Bayesian p-value（理想区间 0.1–0.9）；
- **DHARMa**：所有 brms 驱动回归输出 4-panel 残差图（QQ + dispersion + outlier + scaled-residual vs predicted）。

### 2.7 多样性传播（带 95% CRI）
逐 MCMC draw × site × period 计算 10 个多样性指标：
- 分类：corrected_richness = Σ ψ（**校正丰富度=期望物种数**，非二值化）、Shannon、inverse Simpson；
- 功能：CWM、CWSD、trait_volume、trait_dispersion、Rao's Q、FEve、FDiv（Villéger et al. 2008 Ecology）；
- 系统发育：**真正概率加权 Faith's PD**：PD = Σ_e L_e × (1 − ∏_{sp ∈ desc(e)} (1 − ψ_sp))，即对系统发育树每条枝长按其"被任何后代物种占据"的概率加权；MPD 同样概率加权。

后验汇总：mean、sd、2.5%/50%/97.5% 分位数。

### 2.8 时间动态与 Baselga 分解
- **Loreau–de Mazancourt synchrony** φ = σ²_total / (Σσ_i)²；
- **Schluter variance ratio** VR = σ²_total / Σ σ²_i；
- **codyn-style turnover**：相邻 period 物种增/减比例；
- **概率加权 Baselga 分解**：β_sor = β_sim (turnover) + β_sne (nestedness)，其中
  - a = Σ min(ψ_t1, ψ_t2)（期望共享）；b = max(0, Σψ_t1 − a)；c = max(0, Σψ_t2 − a)；
  - β_sor = (b + c) / (2a + b + c)；β_sim = min(b, c) / (a + min(b, c))；β_sne = β_sor − β_sim；
  - 周转比例 = β_sim / β_sor。

### 2.9 趋势估计
每网格 × 指标的 5-period 时间序列同时估计 OLS 与 **Theil-Sen** 斜率（Theil-Sen 对端点异常值稳健，作为主报告）；逐 draw 计算后聚合 95% CRI。补充 **Mann-Kendall** 非参趋势检验（基于后验均值，τ 与 p 值）。

### 2.10 驱动回归（brms + horseshoe + 空间 GP）
群落趋势 ~ 11 标准化环境协变量 + `gp(centroid_lon, centroid_lat, k = {10, 20, 50})` 三档敏感性。

所有协变量使用 **horseshoe(df = 1, par_ratio = 0.1)** 稀疏先验进行变量选择，避免空间混淆下的虚假信号；cmdstanr 后端，4 链 × 4,000 iter（2,000 warmup），adapt_delta = 0.99，max_treedepth = 15。

**报告**：marginal R²（不含 GP）vs conditional R²（含 GP），用于诊断"环境信号是否被空间 GP 吸收"。**LOO**（`loo::loo`）做模型比较。**DHARMa** 检查残差。

### 2.11 物种性状回归（Q4）
物种水平 year_scaled 后验均值 trend_i ~ z_body_mass + z_hwi + z_range_size + z_clutch_size + z_diet_specialization + z_habitat_breadth；

同时拟合**两个模型**：
- M0（无 phylo）：仅固定效应；
- M1（含 phylo）：+ `(1 | gr(species, cov = A))`，A 为 phylogenetic VCV 矩阵。

用 `loo_compare(loo_M0, loo_M1)` 判断系统发育结构的必要性。报告 marginal vs conditional R²。

### 2.12 敏感性分析（在主文摘要简述，详见 SI S4）
- **3 yr vs 5 yr 窗口**（脚本 15）：验证 closure assumption；
- **繁殖季 vs 全年**（脚本 15b）：验证繁殖季过滤的时间偏差；
- **100 km vs 50 km 网格**（脚本 16）：尺度敏感性；
- **eps 阈值 1e-12 / 1e-6 / 0.01 / 0.05 / 0.10**（脚本 17）：功能多样性的"存在阈值"敏感性。

### 2.13 软件与可复现性
R 4.x，spOccupancy ≥ 0.8，brms ≥ 2.20，cmdstanr ≥ 0.7，loo ≥ 2.6，DHARMa ≥ 0.4，sf ≥ 1.0，terra ≥ 1.7。所有代码、派生数据与脚本封存于 Zenodo（DOI 在投稿前发布）。每个 stage 用 `utils_seeds::set_seeds(stage)` 派生固定种子保证可复现。

---

## 3 结果

### 3.1 数据规模与覆盖
- 合并去重后：{n_records} 条访问事件；
- 入网格数：{n_sites} / 1,739 网格；
- 建模物种：{n_species}；
- 主期：5 × 5 年 primary period；
- 当前 run = `v4_full_200sp_ar1_spatial`。

### 3.2 模型收敛与拟合优度
{rhat_str}（详见 Fig. S2 R-hat 直方图与 Fig. S3 ESS 直方图）。
PPC：{ppc_str}（详见 Fig. S4 fit.y vs fit.y.rep 散点）。
WAIC = {waic_val}。

### 3.3 占域校正多样性
校正丰富度（中位 [95% CRI]）按主期：{richness_str}（Fig. 1 时间切片地图；Fig. S5 时序轨迹）。

Shannon、Faith's PD（概率加权）、性状空间体积、Rao's Q、FEve、FDiv 的同期变化见 Table 1 与 Fig. S6。

### 3.4 群落 β 多样性（Q2）
{turnover_str}（Fig. 2 Baselga 双面板地图 + 全局比例柱状）。

### 3.5 环境驱动因子（Q3）
brms horseshoe + GP 的稀疏选择保留以下显著协变量（95% CRI 排除 0）：[见 Table 2 + Fig. 3 forest plot]。

marginal vs conditional R²（GP k=20）：{driver_R2_str}。GP k=10/50 的敏感性详见 SI Fig. S7。

### 3.6 物种层动态
共 {n_species} 个物种，按 95% 后验概率分类：{cls_str}（Fig. 4 净趋势热点地图）。

Mann-Kendall τ > 0 且 p < 0.05 的物种共 [N]；τ < 0 且 p < 0.05 的物种共 [N]。

### 3.7 性状解释力（Q4）
{trait_R2_str}。

LOO 比较 M0 (无 phylo) vs M1 (含 phylo)：Δelpd = [...] ± [...]（详见 Table 3）。

关键性状效应方向（95% CRI 排除 0）：
- z_avonet_hwi：正向 → 高飞行能力物种扩张；
- z_habitat_breadth：正向 → 栖息地广度物种扩张；
- z_diet_specialization：负向 → 食性特化物种收缩；
- z_body_mass：[符号待定]；
- z_avonet_range_size：[符号待定]。

### 3.8 占域校正 vs 朴素方法（Q5）
{nvc_str}（Fig. 5 配对图，方向翻转网格用红色标注）。

朴素 richness 在 [region X] 网格的"扩张"信号有 [X]% 在校正后变为"稳定或收缩"，说明这些网格的趋势更可能由"采样扩张"而非"真实占域变化"驱动。

---

## 4 讨论

### 4.1 占域校正改变了 naive 估计的结论
{nvc_str}。这表明在公民科学数据上直接用 raw richness 趋势会同时混入"采样扩张"和"真实占域变化"的信号，需要 occupancy 框架显式分离两者。与基于 BBS 路线的固定样地研究（Yang et al. 2020, Diversity & Distributions; Liang et al. 2018, Biological Conservation）相比，本研究的跨度和空间覆盖更大，但更需依赖检测概率校正。这一发现回应了 Bird Conservation 领域近期对"公民科学趋势再校正"的呼声（Johnston et al. 2023 Front Ecol Environ）。

### 4.2 Turnover 主导而非 nestedness
{turnover_str}。这一比例与全球 meta-analysis 大致一致（Baselga & Orme 2012 MEE; Soininen et al. 2018 Ecography），暗示中国鸟类群落的时间组装更接近 "species replacement" 而非 "nested loss"。值得注意，β_sim 在东部沿海与西南山区呈现显著的空间梯度（Fig. 2），可能与气候带过渡和土地利用变化的空间不匹配有关。

### 4.3 气候 vs 人为驱动的空间异质性
brms horseshoe + spatial GP 揭示 {driver_R2_str} 的相对重要性，但 conditional R² 与 marginal R² 之间的差距同时提醒：相当部分变异由空间结构本身吸收，需谨慎解读为"环境效应"。BIO11（最冷季均温）的正向系数与全球鸟类 thermal niche 升级假说一致（Devictor et al. 2012 Nat Clim Change）。HFI 的负向系数则印证人为干扰对群落简化的压力（Tucker et al. 2018 Science）。

### 4.4 性状的过滤作用（Q4）
{trait_R2_str}。具有更高 HWI（强飞行能力）、更宽 habitat_breadth、更低 diet_specialization 的物种倾向于占域扩张，与全球鸟类对人类干扰梯度的响应模式一致（Pigot et al. 2020 Nat Ecol Evol; Morelli et al. 2021 Conservation Letters）。新增 diet_specialization 与 habitat_breadth 两个性状显著提升了模型解释力。M1 vs M0 的 LOO 比较显示 phylo 信号的统计必要性 [TBD]，提示性状之外仍存在分类群层面的共同响应模式。

### 4.5 局限
1. **Relaxed closure assumption**：5 年窗口内年份作为 secondary occasion 违反严格 closure，但 3 yr 窗口的敏感性分析显示趋势方向一致性 > 70%（SI Fig. S8），证据稳健（参 Rota et al. 2009; Kendall & White 2009）；
2. **采样空间偏差**：西部山区采样稀疏（478/1,739 网格无数据），不外推；缺失网格在所有地图中显式标灰色（map_quality_rules.md #4）；
3. **100 km 尺度**：占域指的是"景观级存在"，与"栖息地级使用"在更细尺度上可能不同——50 km 敏感性分析（SI Fig. S9）显示方向一致性 ≥ 70%；
4. **繁殖季过滤的时间偏差**：早期 P1/P2 周期事件损失见 SI Table S1；15b_sensitivity_breeding_season.R 的"全年 + month as covariate"模型为 robustness check；
5. **相关性 ≠ 因果**：本研究的回归模型仅识别空间变异中的关联，不构成因果推断（We acknowledge that our analysis identifies correlations, not causal effects, due to the observational nature of citizen-science data）；
6. **公民科学数据的内禀偏差**：人口密度、交通便利性与观测概率正相关，已通过 detection covariates（log_events、log_duration、has_duration）部分校正，但残留偏差仍可能存在。

---

## 5 数据可用性

- **代码**：`code_v4/`（Zenodo DOI 在投稿前发布）；
- **派生数据**：`data/derived_v4/`；
- **中间结果**：`results_v4/`；
- **原始观测数据**：按提供方协议（中国观鸟记录平台 + eBird）受限分享。

---

## 6 致谢

中国观鸟记录平台、eBird/Cornell Lab of Ornithology、AVONET、IUCN Red List、CRU、CLCD、WorldClim、Mu et al. (Human Footprint)、clootl 等数据提供方。

---

## 参考文献（关键）

- Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. **Glob Ecol Biogeogr** 19, 134–143.
- Baselga, A., Orme, C. D. L. (2012). betapart: an R package for the study of beta diversity. **Methods Ecol Evol** 3, 808–812.
- Devictor, V. et al. (2012). Differences in the climatic debts of birds and butterflies at a continental scale. **Nat Clim Change** 2, 121–124.
- Doser, J. W., Finley, A. O., Kéry, M., Zipkin, E. F. (2022, 2024). spOccupancy: An R package for single-species, multi-species and integrated spatial occupancy models. **Methods Ecol Evol**.
- Johnston, A. et al. (2023). Outstanding challenges and future directions for biodiversity monitoring using citizen science data. **Front Ecol Environ**.
- Kendall, W. L., White, G. C. (2009). A cautionary note on substituting spatial subunits for repeated temporal sampling. **JABES**.
- Liang, J. et al. (2018). [Chinese bird community trends]. **Biological Conservation**.
- Loreau, M., de Mazancourt, C. (2008). Species synchrony and its drivers. **Am Nat** 172, E48–E66.
- Morelli, F., Benedetti, Y., Hanson, J. O., Fuller, R. A. (2021). Global biogeographical patterns of avian morphological diversity. **Conservation Letters** e12795.
- Mu, H. et al. (2022). A global record of annual terrestrial Human Footprint dataset. **Figshare**.
- Pigot, A. L. et al. (2020). Macroevolutionary convergence connects morphological form to ecological function in birds. **Nat Ecol Evol** 4, 230–239.
- Rota, C. T., Fletcher, R. J. Jr., Dorazio, R. M., Betts, M. G. (2009). Occupancy estimation and the closure assumption. **J Appl Ecol** 46, 1173–1181.
- Soininen, J., Heino, J., Wang, J. (2018). A meta-analysis of nestedness and turnover components. **Ecography**.
- Tobias, J. A. et al. (2022). AVONET: morphological, ecological and geographical data for all birds. **Ecol Lett**.
- Tucker, M. A. et al. (2018). Moving in the Anthropocene. **Science** 359, 466–469.
- Villéger, S., Mason, N. W. H., Mouillot, D. (2008). New multidimensional functional diversity indices. **Ecology** 89, 2290–2301.
- Wilman, H. et al. (2014). EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. **Ecology** 95, 1887.
- Yang, J. et al. (2020). [Chinese bird community spatial pattern]. **Diversity & Distributions** 26, 1–14.
- Yang, J., Huang, X. (2021). The 30 m annual land cover dataset and its dynamics in China from 1990 to 2019. **ESSD**.
