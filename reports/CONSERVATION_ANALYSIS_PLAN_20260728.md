# 机制归因 + 保护成效/空缺/规划 —— 分析方案与运行指南

> 目标：把当前论文从"记录了同质化现象"推进到"解释了为什么 + 指出该怎么办"。
> 新增三个脚本（30/31/32），依次回答：**物种数为何增长 → 保护区起没起作用 → 该保护哪里**。
> 日期：2026-07-28

---

## 零、一个先于分析的数学事实

每格物种数是各物种占域概率之和，故全国均值满足恒等式：

**R̄(t) = (S/N) × 平均物种期望分布面积(t)**

S=200、N=1247 固定 ⇒ **每格物种数 +27% 等价于"平均物种分布面积 +27%"**。代入实测：平均占域概率 ψ̄ 从 **0.40 → 0.50**。

所以"物种数增长是否由分布范围扩张引起"——**定义上就是**。真问题是：**哪一种范围扩张、由什么驱动**。脚本 30 的第 2 步会用实际数据自动校验此恒等式。

---

## 一、脚本 30：物种数增长的机制归因

`code_v3/30_range_expansion_mechanism.R`

| 分析 | 输出表 | 判读 |
|---|---|---|
| **扩张 vs 填充分解** | `table_range_decomposition_summary_*` | `pct_expansion` 占优 ⇒ 边缘扩张（气候位移特征）；`pct_infilling` 占优 ⇒ 范围内填充（种群恢复特征） |
| 定殖-灭绝分解 | `table_colonization_extinction_*` | 净增长来自"定殖多"还是"灭绝少" |
| **重心纬度/海拔位移** | `table_centroid_shift_*` | 系统性北移+上移 ⇒ 气候变暖的方向性指纹 |
| **群落温度指数 CTI** | `table_cti_*` | CTI 持续上升 ⇒ 热适应种占比上升（气候驱动经典证据） |
| **生境类群归因** | `table_guild_contribution_*` | 增量若高度集中于少数类群（如水鸟）⇒ **直接解释功能多样性为何没跟上** |

阈值敏感性：核心区阈值取 0.5，另跑 0.3/0.7。

---

## 二、脚本 31：保护成效与保护空缺

`code_v3/31_protected_area_effectiveness.R`

**数据**：已下载 Zenodo《中国自然保护区名录与矢量边界》（CC-BY-4.0，DOI 10.5281/zenodo.14875797），
存于 `data/external/protected_areas/`。实测：**1028 个多边形**，EPSG:4326；
2000 年前建立 **736** 个、2000–2024 建立 **285** 个；类型含森林生态 407、野生动物 319、内陆湿地 109、海洋海岸 33；
级别含国家级 544；总面积 137.4 万 km²（约占国土 14.3%，与公开统计吻合）。

| 分析 | 方法 | 输出 |
|---|---|---|
| 网格覆盖率 | 等积投影下逐格求交，按级别/类型/年代分层 | `table_pa_grid_coverage_*` |
| **保护成效（核心）** | 倾向性得分匹配（海拔、地形起伏、基线 HFI、气候、地被、经纬度），比较区内外多维趋势 ATT | `table_pa_matching_att_*` |
| **类型分解** | 湿地类 vs 森林类保护区分别匹配 | `table_pa_type_effects_*` |
| **双重差分** | 2000–2024 新建保护区（建立前后 × 处理对照） | `table_pa_did_*` |
| **保护空缺** | 物种占域加权保护比例 vs 30% GBF 目标；栖息地宽度与保护比例的相关 | `table_conservation_gap_*` |

**关键判读**：
- `trend_fric_prob` / `trend_rao_q` 的 **ATT 显著为正** ⇒ 保护区延缓了功能同质化（保护成效的核心证据）。
- 若**湿地类保护区效应 > 森林类** ⇒ 支持"湿地恢复驱动水鸟回归"这一机制假说。
- 若"栖息地宽度 vs 保护比例"**rho 显著为正** ⇒ 泛化种反而受更好保护，**特化种存在系统性保护空缺**。

⚠️ 因果强度声明：匹配无法排除未观测混杂，结论应表述为"与保护相关的差异"；DiD 提供更强的准实验证据。矢量仅覆盖 1028/3376 条（偏国家级与大面积），须在文中声明。

---

## 三、脚本 32：系统保护规划与设计

`code_v3/32_conservation_prioritization.R`

在 30% 陆域预算（对标昆明-蒙特利尔 GBF Target 3，另跑 17% 敏感性）下求解三套优先区方案：

| 方案 | 目标函数 | 代表的规划范式 |
|---|---|---|
| **A 物种数导向** | 各物种等权 | 目前主流做法 |
| **B 功能独特性导向** | 按性状空间平均距离加权（越独特权重越高） | 防止功能同质化 |
| **C 同质化风险导向** | 按功能体积下降速率加权网格 | 抢救正在退化的地方 |

**核心产出与政策含义**：
- `table_priority_overlap_*`：三方案两两 **Jaccard 重叠度**。**重叠度越低 ⇒ 只按物种数规划会系统性错过阻止同质化的关键区域**——这可能是全篇最具政策冲击力的结果。
- `table_priority_pa_shortfall_*`：现有保护区对各方案优先区的达成率；`n_priority_unprotected` 即**扩建/新建的首选目标网格**。

---

## 四、运行顺序（23 服务器）

```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3
export V3_CODE_DIR="$PWD"; export V3_RUN_LABEL="v3_full_200sp_ar1_spatial"

Rscript 30_range_expansion_mechanism.R      # 机制归因（需 4D psi）
Rscript 31_protected_area_effectiveness.R   # 保护成效与空缺（需 30 的 CTI 可选）
Rscript 32_conservation_prioritization.R    # 保护规划（需 31 的覆盖率表）
```

依赖：`MatchIt`、`prioritizr` + 求解器（`highs`）、`sf`。三者均需 4D `psi_samples_thinned`。
保护区数据需从本机同步到服务器：`data/external/protected_areas/`。

---

## 五、这些结果将补进论文的哪里

| 结果 | 补入位置 | 作用 |
|---|---|---|
| 扩张/填充分解、CTI、重心位移 | 新增 §3.x 结果 + §4.5 驱动分离 | 把"驱动未定"升级为**有证据的机制归因** |
| 类群归因 | §4.4 机制 | 解释**功能维度为何没跟上**（增量集中于少数类群） |
| 保护成效 ATT / DiD | 新增 §3.x + §4.7 | 直面"中国保护投入"这一竞争性解释（Ouyang 2016；Bryan 2018） |
| 保护空缺 | §4.7 保护启示 | 特化种空缺 = 具体可行动的结论 |
| **规划方案错配** | 新增讨论小节 | **"按物种数规划会错过关键区域"**——从科学发现走向政策建议 |

需新增引用（均已 Crossref 核验 DOI）：Andam 2008 PNAS `10.1073/pnas.0800437105`；Joppa & Pfaff 2009 PLoS ONE `10.1371/journal.pone.0008273`；Ferraro & Hanauer 2014 `10.1146/annurev-environ-101813-013230`；Geldmann 2019 PNAS `10.1073/pnas.1908221116`；Watson 2014 Nature `10.1038/nature13947`；Maxwell 2020 Nature `10.1038/s41586-020-2773-z`；Rodrigues 2004 Nature `10.1038/nature02422`；Margules & Pressey 2000 Nature `10.1038/35012251`；Pollock 2017 Nature `10.1038/nature22368`；Brum 2017 PNAS `10.1073/pnas.1706461114`；Jung 2021 NEE `10.1038/s41559-021-01528-7`；Butchart 2015 Conserv Lett `10.1111/conl.12158`；Xu 2017 PNAS `10.1073/pnas.1620503114`；Hanson 2024 Conserv Biol `10.1111/cobi.14376`；数据集 Zenodo `10.5281/zenodo.14875797`。
