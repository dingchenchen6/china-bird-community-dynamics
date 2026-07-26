# 中国鸟类动态占域研究 — 生态学专家 & 高级编辑独立审查报告

> 审查人：生态学资深审稿人 / 期刊高级编辑视角
> 审查对象：v4 代码、手稿骨架、方法学设计、生态学解释
> 日期：2026-05-11
> 独立性声明：本报告与 REVIEW_v4_full_audit_cn.md（代码自审）互补，侧重方法学合理性与生态学叙事深度

---

## 0 总体评估

| 维度 | 评分 | 说明 |
|---|---|---|
| 科学立意 | A | 全国尺度、25年、stMsPGOcc + 概率加权PD，国内首创 |
| 方法学创新 | A− | 概率加权Baselga + 概率加权Faith's PD 是真创新；stMsPGOcc 本身是工具应用 |
| 数据基础 | B+ | 748万去重事件，国内最大；但系统性偏差（双源、eBird覆盖偏差）未充分处理 |
| 生态学解释深度 | C+ | 核心发现（泛化种扩张）合理但讨论偏薄，缺保护政策启示 |
| 统计严谨性 | B | 大模型设定合理但若干参数可识别性存疑（AR1@5期、NNGP 5邻居） |
| 投稿成熟度 | C+（60/100） | 与自评一致；FULL尚未跑完，数字与最新模型脱节 |

**综合判断**：这是一篇有潜力冲击 Global Change Biology / Nature Communications 的工作，但当前完成度约 60%。按已有路线图推进 6–8 周，配合本报告的补充建议，可达到 85/100 的可投状态。

---

## 1 方法学层面的深层问题（审稿人必问）

### 1.1 AR(1) 时间相关在 5 个主期下的可识别性 [P0]

**问题**：stMsPGOcc 的 `ar1 = TRUE` 在每个站点上拟合 ρ（一阶自相关）和 σ²_t（时间方差）。但仅有 5 个主期 → 相邻时间点仅 4 个间隔。AR1 参数（尤其是 ρ）的估计自由度极低，后验可能非常宽。

**审稿人会问**："With only 5 primary periods, how identifiable is the AR1 temporal correlation? Please show the posterior distribution of ρ and compare models with and without AR1 using WAIC/LOO."

**建议**：
1. 在 `04b_recover_diagnostics.R` 中输出 `rho` 的后验分布（按物种和群落层）
2. 若 ρ 的后验 95% CRI 几乎覆盖整个 (−1, 1) 区间，则 AR1 不可识别，应改用无 AR1 的模型或更简单的随机游走
3. 在 SI 中跑一版 `ar1 = FALSE` 的对比，报告 WAIC 差异
4. 手稿 Methods 中需明确说明："AR1 was included despite limited temporal resolution; posterior checks confirmed identifiable ρ for X% of species"（或坦承不可识别而改用替代方案）

### 1.2 NNGP 5 邻居在中国大尺度是否充足 [P1]

**问题**：100 km 网格、中国东西跨度 ~5000 km，空间相关距离可能达到数百至上千公里。NNGP 仅 5 邻居意味着每个网格只"看到"周边 500 km 内的 5 个最近邻。对于大尺度梯度（如纬度梯度、季风影响），这可能导致空间结构被低估。

**文献参考**：Doser et al. (2022 spOccupancy 论文) 建议 n.neighbors 通常取 5–15，大尺度研究建议 10+。Finley et al. (2019) 的 NNGP 原论文也指出邻居数影响近似精度。

**审稿人会问**："Why 5 neighbors? How sensitive are your spatial occupancy estimates to n.neighbors = 10 or 15?"

**建议**：
1. 在 pilot 中对比 n.neighbors = 5, 10, 15 对 WAIC 和 φ（空间衰减参数）后验的影响
2. 若 WAIC 差异 < 2，则 5 邻居可接受；否则需在 Methods 中论证
3. 检查 φ 的后验：若 φ 集中在先验下限 0.1（对应有效距离 ~10 km，远小于 100 km 网格），说明空间效应几乎不可识别——这会是致命问题

### 1.3 phi 先验 Unif(0.1, 10) 的尺度合理性 [P1]

**问题**：phi 的单位是 km⁻¹。Unif(0.1, 10) 对应的有效空间相关距离范围是 0.1–10 km。但网格分辨率是 100 km！如果真实的空间相关距离是 300–500 km（对应 phi ~ 0.002–0.003），那 phi 的先验完全排除了真实值。

**验证**：在 `04_run_stMsPGOcc_main.R` 中，phi_init = 3 / mean(pairwise_distance)。若平均 pairwise distance ~2000 km，则 phi_init ~ 0.0015，远低于先验下限 0.1。这意味着初始化值就在先验之外！MCMC 会被强行拉回 [0.1, 10] 区间。

**后果**：空间效应被严重低估，NNGP 退化为近乎独立的站点效应，spatial random effect 几乎失去意义。

**建议**：
1. **立即修正** phi 先验为 `Unif(0.001, 1)` 或根据经验设定（如 `Unif(3/max_dist, 3/min_dist)`）
2. 或改用 "range" 参数化（如 spOccupancy 支持的 range = 3/phi），令先验在距离尺度上均匀
3. 跑完后检查 phi 后验的位置——如果集中在 0.1 附近，说明先验在限制推断，需要重新跑

### 1.4 检测协变量缺少 "source" 效应 [P0]

**问题**：eBird 和中国观鸟记录平台的观测协议、物种鉴定能力、空间覆盖模式存在系统性差异。eBird 以 checklist 为主，有严格的协议规范；中国观鸟平台记录方式更灵活。两源的检测概率可能不同。

**审稿人会问**："Did you account for source-specific detection heterogeneity? eBird and local platforms likely differ in protocols and observer expertise."

**建议**：
1. 将 `source`（eBird vs 国内平台）作为 detection 协变量加入模型：
   `det.formula = ~ log_events + log_duration + has_duration + source`
2. 如果 source 与 grid 存在共线（某些网格只有单一来源），则改用 source-specific effort index
3. 在 SI 中报告 source 效应的方向和大小

### 1.5 n.factors = 5 对 200 物种是否充足 [P1]

**问题**：latent factor 模型的 n.factors 决定种间相关的维度。Kéry & Royle (Occupancy Estimation and Modeling) 和 Doser et al. (2022) 建议使用 sqrt(n_species) 或更多。对于 200 物种，sqrt(200) ≈ 14，当前仅用 5 个因子可能无法充分捕捉种间相关结构。

**后果**：共享的种间响应被过度压缩，可能导致占域估计的物种特异性偏差。

**建议**：
1. 在 pilot 中对比 n.factors = 5, 10, 15 的 WAIC
2. 若 WAIC 下降显著，FULL 跑 n.factors = 10（内存允许的话）
3. 在 Methods 中明确报告选择的 n.factors 及理由

---

## 2 生态学数据层面的根本性问题

### 2.1 双源数据的系统性空间偏差 [P1]

**问题**：eBird 数据高度偏向旅游热点（如云南、四川、北京、上海、广东）和交通便利区域。中国观鸟平台的数据虽然更分散，但早期（2000–2010）同样集中在人口密集区。两种源的合并不能消除这种系统性偏差——它只能被检测协变量部分解释。

**更深层的问题**：如果检测概率的空间变化与真实占域的空间变化存在共线（如城市周边调查更多、城市周边鸟类群落确实不同），占域模型也可能无法完全分离两者。

**建议**：
1. 在 Discussion 的 Limitations 中明确讨论这一点
2. 在 SI 中出图：`fig_source_spatial_bias.png`——按网格展示 eBird vs 国内平台的记录比例空间分布
3. 考虑在 driver 回归中加入 "eBird proportion" 作为协变量，检查其是否显著（若显著则说明源偏差未被完全校正）

### 2.2 居留型（留鸟/候鸟）分组缺失 [P1]

**问题**：繁殖季过滤（4–8月）后，冬候鸟和旅鸟在该时段的记录极少。200 个建模物种中可能包含大量迁徙种，它们的占域动态与留鸟截然不同（如冬候鸟的"收缩"可能只是迁飞路线变化）。

**审稿人会问**："How do migratory status and breeding/wintering ranges affect your occupancy trends? Did you stratify analyses by migratory behaviour?"

**建议**：
1. 从 AVONET 或 IUCN 获取居留型数据（resident / migrant / partial migrant）
2. 在 postprocessing 中按居留型分组报告趋势：
   - 留鸟：占域变化反映真实栖息地变化
   - 夏候鸟：繁殖地变化
   - 冬候鸟：越冬地变化（但繁殖季过滤后样本极少）
3. 在手稿 Results 中增加一段："When stratified by migratory status, resident species showed X% expansion vs Y% for migrants..."

### 2.3 繁殖季过滤对早期数据的影响被低估 [P0]

**问题**：v4 手稿提到 "早期 P1/P2 周期事件损失见 SI Table S1"，但这不是一个可以敷衍过去的问题。如果 P1（2000–2004）的繁殖季记录只占全年记录的 30%，而 P5（2020–2024）占 70%，那"占域增加"可能部分反映的是"繁殖季调查比例增加"而非真实变化。

**建议**：
1. 在 `02_build_survey_history.R` 中按 period × source 输出繁殖季过滤前后的比例
2. 若 P1 保留率 < 50%，在手稿中必须明确说明并讨论其对趋势的潜在影响
3. 在敏感性分析（15b）中不仅跑"全年 + month covariate"，还要跑"仅保留繁殖季但加 period × source 交互"的检测模型

### 2.4 100 km 尺度的生态学含义 [P1]

**问题**：100 km 网格约等于 1° 纬度。在这个尺度上：
- 单个网格可横跨多个生态区（如东部平原一个网格可能包含农田+湿地+城市）
- 对于分布范围 < 100 km 的特有种，占域信号被严重稀释
- 5 年主期内 "closure" 假设几乎必然被违反（物种在 5 年内完全可能 colonize 或 go locally extinct）

**审稿人会问**："At 100-km resolution, what does 'occupancy' ecologically mean? How does this coarse grain affect interpretation for range-restricted species?"

**建议**：
1. 在 Discussion 中增加一段明确说明：
   "At 100-km resolution, our occupancy estimates reflect landscape-level presence probability rather than habitat-level use. Species with ranges smaller than the grid size are effectively undetectable at this scale, and our inference should be interpreted as broad-scale range dynamics rather than fine-scale habitat occupancy."
2. 筛选建模物种时增加 "range_size > 100 km" 的条件（或至少报告有多少物种的分布范围 < 100 km）
3. 在 SI 中讨论 closure assumption 的 relaxed 解释（Rota et al. 2009; Kendall & White 2009）

---

## 3 统计与后处理层面的问题

### 3.1 概率加权 Baselga 的数学性质 [P1]

**问题**：`div_baselga()` 函数用概率加权重新定义了 a, b, c：
- a = Σ min(ψ₁, ψ₂)
- b = max(0, Σψ₁ − a)
- c = max(0, Σψ₂ − a)

这一定义在数学上合理，但需要验证：
1. 是否满足 β_sor = β_sim + β_sne？
2. β_sim 和 β_sne 是否仍在 [0, 1] 区间内？
3. 当 ψ 为概率（非 0/1）时，Baselga 的分解公理是否仍然成立？

**审稿人会问**："Your probability-weighted Baselga decomposition is innovative but requires proof that the partition axioms hold for probabilistic inputs."

**建议**：
1. 在 SI 中增加一个数学附录（S8），证明概率加权 Baselga 满足 Baselga (2010) 的分解公理
2. 或引用已有文献（如果有类似的概率加权 β 分解方法）
3. 用模拟数据验证：对同一对 ψ 向量，计算概率加权 Baselga 和传统 Baselga（二值化后），检查相关性

### 3.2 brms + horseshoe + GP 的空间混淆 [P1]

**问题**：驱动回归中同时放入 11 个环境协变量 + 空间 GP。即使使用 horseshoe 稀疏先验，空间 GP 也可能吸收环境协变量的效应（spatial confounding）。

**审稿人会问**："How do you disentangle environmental effects from spatial structure? The spatial GP may absorb environmental gradients."

**建议**：
1. 在手稿中明确报告 marginal R²（仅固定效应）vs conditional R²（含 GP）
2. 若 conditional R² >> marginal R²（如 0.65 vs 0.15），则需坦承："Much of the spatial variation is captured by the GP rather than the measured environmental covariates, suggesting unmeasured spatial processes or spatially structured confounding."
3. 考虑运行 restricted spatial regression（RSR）或 spatially varying coefficients 作为敏感性分析

### 3.3 性状回归的样本量与功效 [P2]

**问题**：6 个性状预测 200 个物种的占域趋势，n=200, p=6，统计功效尚可。但：
- 性状间可能存在共线性（如 body_mass 与 HWI 相关）
- AVONET 性状有缺失值，missForest 插补引入额外不确定性
- diet_specialization 和 habitat_breadth 的分类精度可能低于连续性状

**建议**：
1. 在 SI 中报告性状相关矩阵（VIF）
2. 对 missForest 插补的不确定性做敏感性分析（如多重插补）
3. 在性状回归中同时报告标准化效应大小（|β| / SD）而非仅方向

---

## 4 生态学解释与叙事深度

### 4.1 核心发现的叙事张力不足

当前的核心发现是：
1. 校正丰富度在增加
2. β 多样性由 turnover 主导
3. 泛化种（高 HWI、宽 habitat、低 diet specialization）在扩张

这一发现与全球"biotic homogenization"（生物同质化）文献高度一致。问题是：**这一发现对中国有什么特殊性？** 如果没有中国特异性的洞见，审稿人可能会认为这只是一个"用中国数据验证已知模式"的研究。

**建议——增加中国特异性叙事**：
1. **青藏高原与西部的独特性**：西部采样稀疏但变化剧烈。讨论青藏铁路、川藏公路等基础设施对高原鸟类的影响
2. **东部沿海的城市化 vs 湿地保护**：长三角、珠三角的城市化与滨海湿地保护的博弈
3. **南岭-秦岭作为气候避难所**：这些山脉在气候变化下的"物种库"角色
4. **与 Yang et al. (2020, Diversity & Distributions) 的系统对比**：他们在更短时段、更小尺度上发现了什么？本研究在 25 年尺度上如何扩展/修正了他们的结论？

### 4.2 保护政策启示缺失 [P1]

**问题**：作为一个全国尺度的生物多样性研究，完全没有讨论对保护政策的启示。这在 GCB / D&D 等期刊是不可接受的。

**建议增加**：
1. **优先保护区域识别**：哪些网格同时显示出高校正丰富度下降 + 高 turnover + 高 human footprint 增加？这些网格应被标记为"高优先级保护区域"
2. **迁徙走廊的重要性**：如果 turnover 主导，说明物种在移动而非局部灭绝。这暗示维持连通性（connectivity）比维持局部栖息地质量更重要
3. **公民科学监测网络的优化**：基于采样偏差分析，建议未来在哪些区域增加调查（如新疆、青海、西藏的空白网格）

### 4.3 与已发表中国鸟类研究的对比不足

**必须对比的研究**（当前 Discussion 仅提及，缺乏深入对比）：
1. **Yang et al. (2020, Diversity & Distributions)**：他们在 2000s 年代用鸟类分布数据做区系分析——本研究的占域校正结论如何与之对比？
2. **Liang et al. (2018, Biological Conservation)**：中国鸟类群落长期变化——本研究的动态占域结论是否与他们的固定样地结果一致？
3. **Hu et al. (相关研究)**：中国鸟类多样性热点——本研究是否支持/修正了他们的热点识别？

**建议**：在 Discussion 中专门设置一个段落，用表格形式对比本研究与上述研究的关键结论。

---

## 5 文稿与呈现层面的高级编辑意见

### 5.1 标题优化

当前手稿标题（中文）："检测偏差校正下中国鸟类群落的时空动态（2000–2024）"

**问题**：太平淡，没有突出核心发现或方法创新。

**建议**（任选其一）：
1. "隐形的扩张：占域校正揭示中国鸟类群落的检测偏差与真实动态（2000–2024）"
2. "全国尺度鸟类群落正在同质化吗？——检测校正后的 25 年占域证据"
3. "公民科学数据的真相与幻象：stMsPGOcc 框架下的中国鸟类群落重构"

英文建议：
1. "Hidden homogenization: detection-corrected occupancy reveals how Chinese bird communities have reshaped (2000–2024)"
2. "Citizen science tells two stories: why naive richness trends mislead and what occupancy correction reveals about China's birds"

### 5.2 摘要的关键改进

英文摘要结构合理，但需要更具体的数字占位。当前所有发现都是 `{richness_str}` 等模板字段，FULL 跑完后需要填入具体数字。

**建议**：
- "Corrected richness increased by X% [95% CRI: Y–Z] from P1 to P5"
- "Turnover accounted for X% of β-diversity change, with nestedness contributing only Y%"
- "Naive richness trends flipped direction in X% of grids after occupancy correction"
- "Species with HWI > median were X times more likely to show expanding occupancy"

### 5.3 Cover Letter 改进

当前 cover letter 质量良好，但需要：
1. **明确 target journal**——不要写 "Global Change Biology / Nature Communications / Methods in Ecology and Evolution"，这会让编辑觉得你没有认真选择
2. **推荐审稿人可以更精准**：
   - Andrés Baselga 很好
   - 建议增加一位**中国鸟类占域/群落动态专家**（如中山大学的刘阳团队、中科院动物研究所的相关研究者）
   - 建议增加一位**公民科学数据偏差校正专家**（如 Alison Johnston）
3. **突出"概率加权 PD"的方法创新**——这是真正的数学贡献，应在 cover letter 中单独强调

### 5.4 图表叙述的改进

当前正文对手稿中 Figure 1–6 的描述非常简略。需要：
1. 每个 Figure 在 Results 中至少有一段专门描述
2. Figure caption 需要自包含："What", "How", "Key takeaway"
3. 确保所有地图都满足已有审查中的 7 条硬规则

---

## 6 期刊投稿策略

### 6.1 目标期刊排序

**首选：Global Change Biology**
- 契合度：气候变化 + 生物多样性 + 人类干扰 → 完美契合 GCB 的 scope
- 竞争力：方法上足够新颖（全国首次 stMsPGOcc），生态学故事需要加强讨论深度
- 风险：如果被认为"描述性"（descriptive）而非"机制性"（mechanistic），可能被拒

**次选：Nature Communications**
- 契合度：方法创新（概率加权 PD + Baselga）+ 大数据 + 大尺度
- 竞争力：需要更强的"普适性"叙事（为什么中国结果对全球有启示）
- 风险：NC 对生态学机制的要求比 GCB 更高，讨论偏薄是硬伤

**备选：Diversity & Distributions**
- 契合度：生物地理学 + 分布动态 → D&D 的核心 scope
- 竞争力：方法学上绰绰有余，讨论深度也足够
- 风险：影响力因子低于前两者

**方法学备选：Methods in Ecology and Evolution**
- 契合度：如果主打"概率加权多样性指标"的方法创新
- 竞争力：方法创新足够强，但生态学发现可能被认为"附带"

**建议**：首投 Global Change Biology。若被拒（大概率因讨论深度不足），修改后转投 Nature Communications（加强全球普适性叙事）。

### 6.2 审稿人可能提出的 10 个关键问题（及预案）

| # | 问题 | 严重性 | 应对 |
|---|---|---|---|
| 1 | "5 个时期拟合 AR1 是否充分？" | 高 | 跑 ar1=FALSE 对比，报告 WAIC；若 ρ 后验过宽则去掉 AR1 |
| 2 | "NNGP 5 邻居在中国大尺度是否充足？" | 高 | 跑 n.neighbors = 10, 15 敏感性；检查 phi 后验 |
| 3 | "phi 先验是否排除了真实的空间尺度？" | 高 | 修正 phi 先验；重跑后检查 phi 后验 |
| 4 | "检测协变量是否充分？" | 中 | 加入 source 作为 detection covariate |
| 5 | "概率加权 Baselga 的数学性质？" | 中 | SI 数学附录证明 |
| 6 | "双源数据的系统性偏差如何处理？" | 中 | 在 Discussion 中明确讨论；SI 出 source bias 图 |
| 7 | "为什么没有按居留型分组？" | 中 | 补充居留型分层分析 |
| 8 | "100 km 尺度的生态学含义？" | 中 | Discussion 增加 "landscape-level occupancy" 解释 |
| 9 | "空间 GP 是否吸收了环境效应？" | 中 | 报告 marginal vs conditional R²；若 GP 主导则坦承 |
| 10 | "这对保护政策有什么启示？" | 高 | Discussion 增加保护政策段落 |

---

## 7 独立优先级建议（与自审报告互补）

### P0（投稿前必须解决）

1. **[phi 先验修正]** 将 `phi.unif = list(a = 0.1, b = 10)` 改为 `Unif(0.001, 1)` 或基于数据驱动设定。当前先验几乎必然导致空间效应被低估。**这是可能让整篇论文结论被推翻的问题。**
2. **[AR1 可识别性检查]** FULL 跑完后检查 rho 后验。若 95% CRI 覆盖 >80% 的 (−1, 1) 区间，则去掉 AR1 重跑。
3. **[source 加入 detection]** 将 source（eBird vs 国内平台）作为 detection 协变量加入模型，验证双源系统性偏差是否影响占域估计。
4. **[概率加权 Baselga 数学附录]** 在 SI 中增加数学证明，说明概率加权 Baselga 分解满足 Baselga (2010) 的公理。

### P1（投稿前强烈建议）

5. **[居留型分层]** 按 resident / migrant / partial migrant 分组报告趋势，避免迁徙种的"假收缩"混淆结论。
6. **[n.factors 敏感性]** 在 pilot 中对比 n.factors = 5, 10, 15，若 WAIC 差异 > 2 则调整。
7. **[保护政策启示]** 在 Discussion 中增加一段，基于占域收缩热点和驱动因子结果，提出具体保护建议。
8. **[与已发表中国鸟类研究系统对比]** 用表格形式对比 Yang et al. (2020), Liang et al. (2018) 等研究的核心结论。
9. **[source 空间偏差图]** 在 SI 中出图展示 eBird vs 国内平台的空间覆盖差异。

### P2（修改阶段可做）

10. **[NNGP 邻居数敏感性]** 跑 n.neighbors = 10, 15 的 pilot，验证空间估计对邻居数的稳健性。
11. **[marginal vs conditional R² 报告]** 在 brms 驱动回归中明确报告两者差异，若 GP 主导则坦诚讨论空间混淆。
12. **[IUCN 受威胁等级纳入]** 将受威胁等级（CR/EN/VU/NT/LC）作为额外协变量纳入性状回归。
13. **[标题优化]** 改用更有吸引力的标题，突出核心发现（如"biotic homogenization"）。

---

## 8 验证清单（投稿前逐项过）

在已有自审 8 条验证清单基础上，补充以下生态学/方法学验证：

### 8.1 空间模型验证
```bash
# 检查 phi 后验是否被先验束缚
Rscript -e '
  library(qs)
  fit <- qread("data/derived_v4/stMsPGOcc_fit_v4_full_200sp_ar1_spatial.qs")
  phi_post <- fit$phi.samples
  print(quantile(phi_post, c(0.025, 0.5, 0.975)))
  # 若 97.5% 分位数接近 0.1（先验下限），说明空间效应不可识别
'
```

### 8.2 AR1 可识别性验证
```bash
# 检查 rho 后验宽度
Rscript -e '
  library(qs)
  fit <- qread("data/derived_v4/stMsPGOcc_fit_v4_full_200sp_ar1_spatial.qs")
  rho_post <- fit$rho.samples
  print(quantile(rho_post, c(0.025, 0.5, 0.975), na.rm = TRUE))
  # 若 95% CRI 宽度 > 1.5，则 AR1 不可识别
'
```

### 8.3 双源偏差验证
```bash
# 检查 source 效应的显著性（加入 detection 后）
Rscript -e '
  # 读取新的 table_model_summary，确认 source 系数 95% CRI 不包含 0
'
```

### 8.4 概率加权 Baselga 验证
```bash
# 模拟验证：对随机 psi 向量，检查 β_sor = β_sim + β_sne
Rscript -e '
  source("code_v4/utils_diversity.R")
  psi1 <- runif(50); psi2 <- runif(50)
  b <- div_baselga(psi1, psi2)
  print(c(b$beta_sor, b$beta_sim + b$beta_sne, b$beta_sor - (b$beta_sim + b$beta_sne)))
  # 第 3 个数应接近 0（考虑数值精度）
'
```

---

## 9 总结

本项目在**科学立意、方法学创新、数据规模**三个维度上都达到了国际一流水准。概率加权 Faith's PD 和概率加权 Baselga 分解是真方法学贡献，全国尺度的 stMsPGOcc 应用在中文文献中尚无先例。

但当前版本存在**三个可能颠覆结论的方法学问题**：
1. **phi 先验排除了真实的空间尺度** → 空间效应可能被严重低估
2. **AR1 在 5 个时期下可能不可识别** → 时间结构可能被错误设定
3. **双源数据的系统性偏差未充分处理** → 检测校正可能不彻底

以及**两个生态学叙事的短板**：
4. **保护政策启示完全缺失** → GCB/NC 级别期刊的硬性要求
5. **与已发表中国鸟类研究缺乏系统对比** → 未能建立学术对话

建议按以下顺序推进：
1. **立即修正 phi 先验并重新跑 pilot**（1–2 天）
2. **修正后跑 FULL**（1–2 周，服务器）
3. **AR1 可识别性检查**（FULL 跑完当天）
4. **source 加入 detection + 居留型分组**（1 周）
5. **Discussion 扩展 + 保护政策段落**（3–5 天）
6. **系统对比已发表研究**（3–5 天）
7. **图表自检 + 投稿准备**（1 周）

**预计总时间**：4–6 周（如果 phi 先验修正后不需要重跑 FULL）。

---

> 本报告独立于 REVIEW_v4_full_audit_cn.md（代码层审查），两者互补使用。代码 bug 修复请参照自审报告；方法学与生态学层面的深层问题请参照本报告。
