# 中国鸟类群落动态研究 —— 从头到尾系统审查与梳理总结

> 审查者：资深生态学专家 + 高级编辑视角
> 日期：2026-06-02
> 审查方法：**回到真实结果表（`results_v3/*.csv`）与核心代码（`code_v3/*.R`）逐环节 ground-truth**，不依赖此前已过时的 `PROJECT_AUDIT_500SP_EVIDENCE`（该审计写于本项目 4 链重跑之前）。
> 结论一句话：**这是一个数据规模大、方法学扎实、结论有冲击力的研究，成熟度远高于早先审计所述——4 链收敛、后验 CrI、混淆控制、敏感性分析绝大多数已完成；主结论 H1/H3 统计稳健；核心"同质化"叙事(H2)方向真实但统计闭环尚差三步。**

---

## 一、执行摘要：研究现在处于什么位置

早先的 `PROJECT_AUDIT`（单链、无诊断、无 fit 对象）**已全面过时**。核对当前 `results_v3/` 真实产物后确认，项目已在服务器完成 4 链重跑与完整后处理，证据链条基本闭合：

| 环节 | 早先审计说 | **真实现状（已核对）** |
|---|---|---|
| MCMC 链数 | 单链 | **4 链**（5000 burn-in, thin 2, 1250/链）`table_*_4chain.csv` |
| 收敛诊断 | 无 | **有** R̂/ESS 表（200sp & 500sp）|
| 后验 CrI | 无 | **有** `table_community_metrics_with_cri_*`（逐格逐期 mean/l95/u95）|
| effort 混淆控制 | 未做 | **有** `table_effort_confound_controls_*`（4 项）|
| 敏感性分析 | 缺 | **有** 3yr 窗口/繁殖季/网格尺度/eps 阈值脚本(15–17)|
| PPC | 缺 | 脚本在(05c)，数值待确认导出 |

**主结论稳健性判断：**
- **H1（占域校正后物种数上升）——稳健。** richness 79.15→100.83（5 期单调），Shannon/invSimpson 同向上升，4 链 CrI 支持。
- **H3（晚期 beta 转向 nestedness）——稳健。** Baselga turnover 65.9%→11.3%，P4–P5 的 95% CrI 上界仅 **0.166**，统计上确为 nestedness 主导。
- **H2（功能不随物种数增长／隐性同质化）——方向真实、统计闭环未完成。** 这是全篇卖点，也是当前最脆弱一环（详见三、④⑤）。
- **H5（检测校正改变推断 + 扩张稳健于观察者扩张）——大部分成立**，但混淆控制仍差 source-term 一项，且"固定-effort 反事实"证据性质需精确表述（已在稿件修正）。

---

## 二、逐环节系统审查

### ① 数据整合 —— 稳健
CBRC + GBIF + eBird，2000–2024，100 km 网格 **1247 格**，5 期。source-aware 去重（species×date×coord×source×observer）。`table_01_merge_summary`、`table_02_survey_summary` 齐备。✔ 无重大问题；投稿需附去重审计与逐源年度覆盖表（稿件已列 S1）。

### ② 占域模型 —— 稳健，一处需报告
- **200sp `stMsPGOcc`**（空间 NNGP + AR1），4 链。psi 后验为 4D（draws×species×sites×periods，05 脚本 FIX #3 已保证）。
- **500sp `tMsPGOcc`**（时间、非空间、10 latent factors），4 链，仅作广度检验（稿件全程标注 non-spatial ✔）。
- ⚠️ **收敛细节**：R̂ 全部 <1.05（max 1.039，为 `tau.sq.beta[11]` 方差超参）；**ESS 最小 70，38 个群体级参数中 3 个 ESS<100、17 个<400**。R̂ 合格；**ESS 偏低集中在社区级方差超参（tau.sq.beta），属常见的慢混合，可容忍但必须如实报告**，并建议补报占域概率/派生量层面的 ESS（真正的推断依据）。稿件原"min bulk ESS=70"表述不准（诊断表为单一 `ess` 列、非 bulk/tail，且仅 38 个群体级参数）——**已在稿件修正**。

### ③ 派生多样性指标 —— 稳健（本地已复核）
从 `table_community_metrics_with_cri` 复算全国各期均值，与稿件数字一致：
- richness 79.15→82.23→87.02→90.53→**100.83**（+5.16/期，单调）
- Shannon 4.74→4.88（+0.034/期）；invSimpson 106.4→126.6
- **trait_volume 1.383→1.378→1.376→1.368→1.365（−0.0047/期，5 期单调降）**
- **rao_q 1.705→1.699→1.696→1.684→1.680（−0.0066/期，单调降）**
- feve +0.0017/期、fdiv +0.0020/期（微升）
→ 功能体积与 Rao's Q 的**下降是单调的、方向真实**，不是随机抖动；但绝对幅度小（−1.3%/−1.5%），**严格的 P(decline) 需从 trend draws（05e 产物）计算**（见三、④）。

### ④ Beta 多样性与同质化 —— 一稳一缺
- **时间 beta（Baselga turnover/nestedness）——稳健、带 CrI。** 见上 H3。✔
- 🔴 **空间 beta（跨格点相似度随时间，同质化"金标准"）——缺失且因代码 bug 全 NA。** `table_homogenization_trend_*` 5 期 `mean_sorensen` 全为 NA。**根因已定位**：`12_homogenization_spatiotemporal_maps.R` 第 92 行只对 4D psi 计算，遇 3D psi 静默 `next` 跳过全部——历史 thinning 版本的 3D psi 正是致 NA 之因。**已修复**（见四）。

### ⑤ effort 混淆控制 —— 大体成立，差一项 + 一处逻辑
真实表 `table_effort_confound_controls` 核对：
- detection~log_events slope **0.925 [0.897–0.953]**：effort 强烈进入检测层 ✔
- **effort 饱和子集**：corrected-richness slope median **4.07**、mean 4.54、[2.36–7.86]、**n=30 格**。CI 排除 0（支持真实扩张），但 **n=30 样本小、代表性弱**，且原稿用 mean 未用更稳健的 median——**已在稿件补 median 并标注"小而 effort 稳定子集"**。
- **固定-effort 反事实**：表中量 = detection_prob_change_P1→P5 = **0.323 [0.304–0.340]**。此数**证明的是"effort 增长使检测概率升 0.32、被检测层吸收"（必要条件）**，而非"固定 effort 后 occupancy 趋势仍在"（充分条件）。原稿逻辑错配——**已在稿件修正为如实表述**。→ 真正击破 confound 的直接证据目前主要靠 n=30 的饱和子集，建议扩展（见下）。
- 🔴 **source-term ΔWAIC = NA（未跑）**：平台异质性检测项还没做（需 model refit）。稿件标 PENDING 正确。

### ⑥ 驱动分析 —— 稳健、定性正确
variance partitioning（baseline 空间/环境 adjR²=0.14 > 气候变化 0.04 > 土地利用 0.01）+ RF importance（经度、基线冬温、极端温变、海拔居前），与 `table_varpart_*`、`table_rf_importance_*` 一致。稿件已恰当地表述为**空间结构化关联而非因果**。✔

### ⑦ 物种趋势与 winner/loser —— 稳健
136 扩张/62 稳定/2 收缩（200sp）；trait/env Spearman（habitat breadth ρ=0.38 等）与 `table_trait_regression_*`、`table_env_trend_correlation_*` 一致。winner=水鸟/开阔地/人伴生，loser=山地林特化种，与"广栖息地泛化种驱动"机制自洽。✔

### ⑧ 敏感性 —— 已具备，建议纳入 SI
3yr 窗口、繁殖季、网格尺度、eps 阈值脚本(15–17)均在。投稿前把结果汇入 SI 即可，无需新做。

---

## 三、分级问题清单

### 🔴 真缺口（需服务器重跑/补算，决定"同质化"主张能否闭环）
1. **空间 beta 全 NA** —— bug 已修（四），需在 server 用 4D psi 重跑 `12` 得到"跨格点相似度随时间下降"这一同质化直接证据。
2. **功能指标 slope 的 P(decline)** —— 从 `05e_export_trend_draws.R` 的 trend draws 计算 trait_volume/rao_q 斜率后验，给出 P(slope<0)；据此锁定标题级 "functional homogenization" 的措辞强度。
3. **source-term ΔWAIC** —— 补 detection 含 source 协变量的 model refit。

### 🟠 严谨性/表述（部分已在本次修正）
4. ESS=70 偏低（社区级方差超参）——如实报告 + 建议补报占域/派生量 ESS；慢混合超参可增迭代。**稿件表述已修正**。
5. 固定-effort 反事实的证据性质——**稿件逻辑已修正**；建议补一个"固定 effort 情景下的 occupancy 趋势预测"真反事实以彻底封堵 Devil's Advocate。
6. effort 饱和子集 n=30——**稿件已补 median 与样本量警示**；若可行，放宽"饱和"定义以增大子集。
7. WAIC 绝对值（984,170）无独立意义——仅在模型比较时报告；**稿件已移除该裸值**。

### 🟢 已稳健（可给作者信心，无需再动）
数据整合、4 链收敛（R̂）、richness/Shannon/turnover 的 CrI、驱动分析、winner/loser、trait/env 关联、敏感性脚本——均已核对、内部自洽。

---

## 四、本次已动手的完善优化

1. **修复 `code_v3/12_homogenization_spatiotemporal_maps.R`（第 92–97 行）**：将 3D psi 的**静默 `next`**改为**显式 `stop` 报错**（含维度与重跑指引），杜绝再产出误导性的全 NA 空间 beta 表；4D psi 路径不变、可正常计算。→ 服务器以当前 4D psi 重跑 `12` 即可得到 spatial-beta 时间序列。
2. **`MANUSCRIPT_GCB_v6` 表述修正两处**：
   - 收敛：`min bulk ESS=70 … WAIC=984,170` → 精确为"max R̂=1.039（38 个群体级参数）、min ESS=70、3 个方差超参 ESS<100"，移除无意义 WAIC 裸值。
   - effort 混淆：把"固定-effort 反事实"数字如实表述为"effort→detection 吸收"，并给 effort 饱和子集补 median(4.07) 与 n=30 警示。

（未改动任何真实数值；未触碰你手动回填的正确数字。）

---

## 五、总体判断与下一步

**判断**：研究本体（设计、数据、模型、主结论）已达顶刊竞争力；H1/H3 稳健，驱动与机制自洽。**唯一挡在"可发表"前的是核心卖点 H2（功能同质化）的统计闭环**——功能下降方向真实（本地已核实单调），但需要①空间 beta（bug 已修，待重跑）与②功能斜率 P(decline)（待 trend draws）两项把"下降"从"点估计单调"升级为"统计显著"。第三项 source-term 是稳健性补强。

**建议下一步（按优先级）**：
1. 服务器重跑 `12`（修复后）→ 得空间 beta；从 trend draws 算 trait_volume/rao_q 的 P(decline)。这两步一完成，"functional homogenization" 即可从 PENDING 措辞升级为定稿强主张。
2. 补 source-term refit（ΔWAIC）。
3. 把敏感性(15–17)结果汇入 SI；补报占域/派生量 ESS。
4. 完成后，我可据实测值把 GCB 稿件所有 `[PENDING]` 精确回填并锁定措辞档位。
