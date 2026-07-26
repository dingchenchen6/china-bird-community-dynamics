# Editorial & Peer-Review Report — Submission Readiness for *Nature Ecology & Evolution*

> Manuscript under review: `MANUSCRIPT_TOP_JOURNAL_REWRITE_20260602.md`
> Title: *Expansion Without Differentiation: Occupancy-Corrected Citizen Science Reveals Functional Homogenization of Chinese Bird Communities, 2000–2024*
> Reviewer role: Acting Senior Editor + 3 subject-expert referees + Devil's Advocate
> Target venue standard: NEE Article (broad significance, ~3,000–3,500-word main text, ≤150-word single-paragraph abstract, online Methods)
> Date: 2026-06-02
> 评审语言：中文说明 + 英文术语；本报告用于驱动后续系统改稿

---

## Part 0 — 编辑总评与决定建议 / Editor's Summary & Decision

**一句话判断：** 这是一篇**科学叙事达到顶刊量级、但当前稿件形态会被 NEE 当场退稿（desk reject）**的工作稿。核心发现（占域校正后物种数上升、但功能体积与 Rao's Q 下降、晚期 beta 转向 nestedness ⇒ "扩张但无功能分化 / 隐性同质化"）是一个真正有冲击力、对 citizen-science 时代生物多样性推断有方法论警示意义的故事。但要进入 NEE 的评审流程，必须先跨过三道门槛。

**若以现状投 NEE，预测结果：Desk reject（不送外审）。** 三个直接触发点：

1. **稿件第 5 行自带 "Evidence status: ... must replace single-chain model summaries with verified multi-chain convergence diagnostics"** —— 等于在投稿信里告诉编辑"我的核心模型还没收敛验证"。NEE 编辑不会送审。
2. **Results 通篇是工作备忘语**（`the submission version should…`、`Current PD panels should be repaired`）—— 这是 lab notebook，不是 manuscript。
3. **Abstract ~300 词 / Introduction 无量化 gap、无假设、无广义意义钩子** —— 不符合 NEE 对"broad audience significance"的硬要求。

**结论：本稿应定位为 "Major revision before it can enter review"。** 潜力评级高，但需一次**实质性重做（模型重跑 + 全文重写）**，而非润色。下面五维打分针对"当前稿件形态"，括号内为"修复后可达"的预估。

### 五维评分 / Five-Dimension Score（满分 5）

| 维度 | 权重 | 当前分 | 修复后潜力 | 关键依据 |
|---|---|---|---|---|
| Originality 原创性 | 20% | 4.0 | 4.5 | 把"detectability 校正"与"多维同质化"在国家尺度连起来，是真正的新意；但需更强的全球文献定位才能撑起 NEE 的"广义重要性" |
| Methodological Rigor 方法严谨性 | 25% | **2.0** | 4.0 | 单链 MCMC、无收敛诊断、无 fit 对象、phi 先验未做尺度校验、检测模型未含 source、500 种被误置于空间叙事旁 —— 当前是**致命短板** |
| Evidence Sufficiency 证据充分性 | 25% | **2.5** | 4.0 | 功能性指标变化幅度极小（如 trait volume 1.383→1.365，约 −1.3%）却未报告不确定性，homogenization 主张目前**统计上脆弱**；PD 主图损坏 |
| Argument Coherence 论证连贯性 | 15% | 3.0 | 4.5 | "扩张但无分化"叙事清晰，但 confound（观察者扩张本身）未被正面击破，因果话术与相关证据之间张力大 |
| Writing Quality 写作质量 | 15% | 2.5 | 4.5 | 结构骨架好，但混入备忘语、Intro/Abstract 未对标 NEE、参考文献仅 9 条 |
| **加权总分** | | **≈ 2.65 / 5** | **≈ 4.25 / 5** | 当前低于送审线；修复后具备 NEE 竞争力 |

---

## Part 1 — Senior Editor 评估（适配性·格式·重要性）

### 1.1 NEE 适配性 / Fit
- **正面：** 主题（生物多样性变化的真实性 vs 观测假象）正是 NEE 关注的"对全领域有方法论与认知意义"的问题；中国 25 年、国家尺度、citizen science，地理与政策语境独特。
- **风险：** NEE 偏好**机制+概念推进**，不只是"我们对一个新地区做了占域模型"。必须把卖点从"我们校正了中国鸟类数据"上移到"**citizen-science 扩张会系统性地把'看见得更多'误读为'生物多样性恢复'，而多维分解揭示这是隐性同质化**"这一普适命题。

### 1.2 格式硬缺口 / Format gaps（NEE Article）
| 项 | NEE 要求 | 当前稿 | 行动 |
|---|---|---|---|
| Abstract | 单段、≤150 词、无引用、无小标题 | ~300 词、2 段 | 重写为单段 ~150 词 |
| Introduction | ~3 段、面向广义读者、结尾给明确问题/假设 | 4 段、缺量化 gap 与假设 | 重构，加 H1–H5 |
| Results 小标题 | 可有，但需"陈述发现"而非"备忘" | 含 `submission version should…` | 全部清除备忘语 |
| Methods | 置于文末/在线，可长 | 位置 OK | 补 source 检测项、phi 先验、概率化 Baselga 推导引到 SI |
| Display items | Article 通常 ≤ ~6 | 7 张主图 plan | 合并到 5–6，PD 图必须修复或换 McTavish |
| References | 顶刊通常 40–70 条 | **仅 9 条** | 扩到 ~45–60，补全球 homogenization / citizen-science / occupancy 文献链 |
| 必备声明 | Data availability, Ethics, Author contrib (CRediT), CoI, Funding, AI disclosure | 仅 Data/Code + Acknowl. | 补齐全部 |

### 1.3 编辑视角的三个 make-or-break
1. **证据可信度（见 Part 2/R1）** —— 没有多链收敛诊断，方法分进不了 4 分区。
2. **同质化主张的统计稳健性（见 Part 3/R2）** —— ~1% 的功能指标变化必须带后验不确定性，否则评审会说"这是噪声"。
3. **观察者扩张 confound 的正面回应（见 Part 5/Devil）** —— 这是全篇的存亡问题，必须有专门段落+敏感性分析正面击破。

---

## Part 2 — Reviewer 1：贝叶斯/占域建模专家

**总体：** 模型选择（stMsPGOcc 空间多物种动态占域 + tMsPGOcc 广度扩展）是恰当且前沿的，但执行与报告达不到顶刊门槛。

**Major:**
1. **(致命) 单链、无收敛诊断。** 两个 summary 表均 `n_chains=1`，无 R-hat / ESS / WAIC / PPC。**必须**重跑 ≥4 链并报告 max R-hat、min bulk/tail ESS、关键参数 trace、posterior predictive check。在收敛证据齐备前，任何点估计都不可作为 NEE 结论。
2. **(致命) 无完整 fit 对象。** 后验无法复现，违反 NEE 的可复现性与数据可得性要求。投稿前必须重生成并归档（Zenodo + checksum）。
3. **500 种模型是时间非空间（tMsPGOcc，无 NNGP）。** 稿件已较谨慎，但叙事上仍与空间结论并置，易误导。必须在每次出现处显式标注"temporal breadth extension, not spatial"，且**不可**用它支持任何空间/地图/空间驱动推断。
4. **phi 空间先验尺度未校验。** `phi.unif=c(0.1,10)` 配 100 km 网格，下界可能过高（对应有效空间程过短），会人为压制空间自相关。需先验敏感性分析 + 报告后验是否撞先验边界。
5. **检测模型缺 source 项。** 仅 `log_events + log_duration + has_duration`；eBird/GBIF 与中国观鸟记录中心协议差异大。若数组支持，**必须**加 source（或 source×effort）检测协变量，否则平台异质性会被吸收进 occupancy —— 直接威胁核心主张（见 Devil）。

**Minor / 需说明:**
6. **闭合假设（relaxed closure）。** 5 年一期、期内年份作重复访问，违背严格闭合；需用更短窗口（如 3 年）做敏感性，并引 Rota et al. 2009 正面讨论。
7. **概率化 Baselga 分解** 用 occupancy 概率向量的 pairwise minima 定义期望共享占域 —— 这是非标准推广，需在 SI 给出有界性与分解性质的数学证明，否则方法学评审会质疑 turnover/nestedness 比例的可解释性。
8. **temporal synchrony = 1.0**（审计指出）疑为实现/尺度 bug，**不要**在正文强调，修复或移除。
9. 趋势分类（expanding/stable/contracting）的阈值与可信区间规则需明确写出；5 期序列建议 Theil–Sen + Mann–Kendall 而非 OLS。

---

## Part 3 — Reviewer 2：宏生态 / 生物多样性变化 / 同质化专家

**总体：** "expansion without differentiation"是漂亮的概念钩子，但目前**证据强度配不上结论强度**。

**Major:**
1. **功能性变化幅度过小且无不确定性。** Trait volume 200sp 1.383→1.365（−1.3%）、Rao's Q 1.705→1.680（−1.5%），500sp 同样 ~1%。**没有后验可信区间，无法判断这是真实下降还是零附近抖动。** homogenization 是全篇的"反转"卖点，必须给出：(a) 每期后验分布/CI；(b) 趋势斜率的后验概率（P(slope<0)）；(c) 与 null/randomization 的对比。否则评审会判定"功能下降"不成立，整个"分化缺失"叙事崩塌。
2. **同质化需要空间维度证据，而非只看 alpha 指标均值。** 真正的 biotic homogenization 通常以**空间 beta 多样性下降 / 群落间相似度上升**来定义。当前主要是时间 beta（相邻期）+ alpha 功能指标。建议补一个**空间 beta 多样性随时间下降**的直接证据（跨格点群落相似度上升），这才是 homogenization 的金标准信号。
3. **turnover→nestedness 的解释有歧义。** 晚期 nestedness 主导既可解释为"广布种到处累积"（你的解读），也可解释为"局域丢失导致嵌套"。需用 winners/losers 的功能身份 + 空间格局区分两者；目前论证停在"consistent with"。
4. **与全球文献的对话几乎为零。** 必须接入：McKinney & Lockwood 类同质化经典、Magurran 时间序列 biodiversity change、Dornelas et al. (turnover without net loss)、Jarzyna & Jetz（多维多样性变化）、Bowler/Callaghan 等 citizen-science 大尺度趋势。当前 9 条引用无法支撑 NEE 的概念定位。

**Minor:**
5. "winners"里大量水鸟/开阔地/人伴生种，与广栖息地解读自洽 —— 这是好证据，应在正文显式做成"功能身份 vs 趋势"的图，而不仅列拉丁名。
6. 保护启示（conservation implication）目前偏弱，NEE 喜欢"so what for the planet/policy"。

---

## Part 4 — Reviewer 3：Citizen science / 检测偏差专家

**总体：** 论文的方法论立意（用占域模型把 effort 从 change 里剥离）正确且重要，但**当前没有充分证明"剥离成功"**。

**Major:**
1. **naive vs corrected 仅 2–3% 格点翻向，被你解读为"correction 主要改幅度"。** 但反方会说：**翻向比例这么低，恰恰说明 occupancy 趋势与 naive 趋势高度同向，无法排除 occupancy 增长本身就携带了未被检测协变量吸收干净的 effort 信号。** 必须更强地证明检测子模型确实吸收了 effort 的时间增长（如：检测概率随 log_events 的后验斜率、不同 effort 水平下的占域估计稳定性）。
2. **effort 的时空格局与所声称的趋势在地理上共线。** 观鸟活动近 25 年在中国东部/城市/沿海爆发式增长，而这些地方正是你报告 occupancy 上升、cropland/built/HFI 正相关的地方。需要一个**"effort 增长 vs occupancy 增长"的空间错位检验**：在 effort 早已饱和的格点，occupancy 是否仍上升？若是，结论稳；若否，则是 detectability 残差。
3. **source 异质性未建模**（与 R1#5 呼应）：不同平台的 detection 截距/斜率不同，未建模会把平台切换误判为占域变化。

**Minor:**
4. 建议补 effort 校正的"安慰剂检验"：对已知稳定/普查充分的常见种，模型是否给出平趋势。

---

## Part 5 — Devil's Advocate：全篇的存亡之问

> **"如果 25 年间观察努力激增、检测概率上升，而你的占域估计也上升 —— 你如何在统计上排除：所谓的'物种扩张'，其实是检测子模型没吸收干净的 effort 残差，被动态占域过程当成了真实的 colonization？"**

这是 NEE 评审一定会问、且能一票否决的问题。当前稿件**没有正面回答**。要存活，改稿必须提供至少 2 条独立防线：

- **防线 A（effort 饱和子集）：** 在观察努力早已饱和、不再随时间增长的格点子集里重估趋势。若 occupancy 仍上升 → 真实信号。
- **防线 B（检测-占域分离诊断）：** 报告检测协变量的后验效应，证明 effort 的时间趋势被检测层显式捕获；并做"固定 effort 反事实预测"。
- **防线 C（source/平台敏感性）：** 加 source 检测项后核心结论不变。

**第二存亡之问（功能同质化）：** ~1% 的功能指标变化在没有 CI 时**不可声称为下降**。若加上 CI 后趋势跨越 0，则全篇标题级主张"functional homogenization"必须降级为"taxonomic expansion without detectable functional expansion"，措辞要随证据强度调整。

---

## Part 6 — 分级问题清单（可操作）/ Prioritized Action List

### 🔴 P0 — 投稿前必须解决（分析层面，需重跑；不解决无法送审）
- **P0-1** ≥4 链重跑 200sp stMsPGOcc 与 500sp tMsPGOcc，报告 R-hat/ESS/WAIC/PPC；归档完整 fit 对象。
- **P0-2** 检测模型加入 source（平台）项，或论证不可识别。
- **P0-3** phi 空间先验尺度校验 + 先验敏感性；确认后验不撞边界。
- **P0-4** 功能多样性各指标（trait volume / Rao's Q / FEve / FDiv）输出**后验 CI 与趋势斜率后验概率**；据此决定 homogenization 措辞强度。
- **P0-5** effort 饱和子集 / 固定-effort 反事实 的 confound 检验（Devil 防线 A/B）。
- **P0-6** 修复 PD 主图：统一用已验证的 McTavish 概率加权 PD/MPD，弃用损坏的 Jetz 灰图。
- **P0-7** 补"空间 beta 多样性随时间下降"的直接同质化证据。

### 🟠 P1 — 进入评审前强烈建议（写作可做 / 部分需分析配合）
- **P1-1** 清除 Results 全部备忘语（`submission version should…` 等）。
- **P1-2** Abstract 重写为 NEE 单段 ~150 词；Introduction 重构为 3 段 + 量化 gap + H1–H5。
- **P1-3** 参考文献扩到 ~45–60，建立与全球 homogenization / citizen-science / occupancy 方法学的对话。
- **P1-4** 500sp 每处显式标注"temporal breadth extension, non-spatial"。
- **P1-5** 闭合假设：3 年窗口敏感性；引 Rota et al. 2009。
- **P1-6** 概率化 Baselga 分解的数学性质证明移入 SI。
- **P1-7** 主图压缩到 5–6 张，做"功能身份 × 趋势"winners/losers 图。

### 🟡 P2 — 提升竞争力 / minor
- **P2-1** Theil–Sen + Mann–Kendall 替代 OLS 趋势。
- **P2-2** 移除/修复 synchrony=1.0。
- **P2-3** 补齐 Ethics / CRediT / CoI / Funding / AI disclosure 声明。
- **P2-4** 50 km 网格 pilot、繁殖季 vs 全年 敏感性放 SI。
- **P2-5** 强化 conservation/policy take-home。

---

## Part 7 — 改稿将怎么做（下一阶段计划，待你确认）

按你选择的"假设将重跑模型、按最终版写"，第二阶段我会产出一份**对标 NEE 的系统重写稿**，具体：

1. **结构对标 NEE Article：** 单段 ~150 词 Abstract → 3 段广义 Introduction（含量化 gap + H1–H5）→ 主题化 Results（清除全部备忘语）→ Discussion（正面回应 confound、给全球意义与保护启示）→ 在线 Methods（补 source 项、phi 先验、概率化 Baselga 引 SI）。
2. **证据诚信处理：** 真实点估计（richness 79.15→100.83 等来自实跑表）保留；**收敛诊断/CI/confound 检验结果用清晰占位符**（如 `[R̂_max = __ ; ESS_bulk,min = __ from 4-chain run]`、`[functional trait-volume slope: posterior mean __, 95% CrI __ to __]`），绝不编造数值 —— 你回填实跑结果即可。
3. **措辞与证据强度挂钩：** homogenization 主张的强弱表述设两个版本开关，依 P0-4 的 CI 结果择一。
4. **配套产出：** 更新版 Figure plan（5–6 张）、~45–60 条参考文献骨架（含 DOI，逐条可核验）、必备声明区块、NEE cover letter、SI 大纲。

> 请确认是否按此进入第二阶段改稿；如需调整目标刊措辞强度或先补某项 P0 分析，请一并告知。
