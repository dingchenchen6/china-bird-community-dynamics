# 服务器运行说明（供 Kimi 或其他助手直接执行）

> 这份文档是**自包含**的：不需要任何先前对话上下文即可照做。
> 目标：在计算服务器上跑完"机制归因 + 保护成效/空缺/规划"三个分析，产出结果表。

---

## 0. 背景（30 秒读懂要做什么）

研究对象是**中国鸟类 2000–2024 年的群落变化**。已完成的部分显示：检测校正后每格物种数上升 27%（79.2→100.8 种），但**功能多样性没有同步上升**，群落趋于同质化。

本次要跑的三个脚本回答两个尚未解决的问题：

1. **物种数为什么增长？**（脚本 30）
2. **保护区起作用了吗、该往哪里扩？**（脚本 31、32）

---

## 1. 服务器与路径

| 项 | 值 |
|---|---|
| 主机 | `162.105.149.23` |
| 用户 | `dingchenchen` |
| 密码 | **由用户单独提供**（未写入本仓库） |
| 项目目录 | `~/Documents/New project/bird_dynamic_occupancy_analysis` |
| 代码目录 | `<项目目录>/code_v3` |
| 结果目录 | `<项目目录>/results_v3` |

登录：

```bash
ssh dingchenchen@162.105.149.23
```

---

## 2. 前置条件（必须先满足，否则脚本会失败）

### 2.1 关键输入：4 维 psi 后验数组

```bash
ls -la ~/Documents/New\ project/bird_dynamic_occupancy_analysis/data/derived_v3/psi_samples_thinned_v3_full_200sp_ar1_spatial.rds
```

必须是 **4 维**（draws × species × sites × periods）。检查：

```bash
Rscript -e 'x <- readRDS("~/Documents/New project/bird_dynamic_occupancy_analysis/data/derived_v3/psi_samples_thinned_v3_full_200sp_ar1_spatial.rds"); print(dim(x$psi_samples_thinned))'
```

- 输出 4 个数字 → ✅ 可以继续
- 输出 3 个数字 → ❌ 需先重跑 `05_postprocess_diversity_extended.R`（它会保存 4 维版本）

### 2.2 保护区数据（脚本 31、32 需要）

需把本地这个目录整个上传到服务器同名位置：

```
data/external/protected_areas/
├── china_nature_reserves.rar
├── READ_ME.md
└── 全国自然保护区名录+矢量边界/
    ├── 保护区.shp / .shx / .dbf / .prj   ← 关键
    └── 全国保护区数据.xlsx
```

从本地上传（在**本地**终端执行）：

```bash
rsync -avz "$HOME/Documents/New project/bird_dynamic_occupancy_analysis/data/external/protected_areas/" \
  dingchenchen@162.105.149.23:"~/Documents/New project/bird_dynamic_occupancy_analysis/data/external/protected_areas/"
```

数据来源：China Nature Reserve Specimen Resource Sharing Platform (2024). *List and Vector Boundaries of Nature Reserves in China*. Zenodo. https://doi.org/10.5281/zenodo.14875797 （CC-BY-4.0）

### 2.3 R 包

```bash
Rscript -e 'install.packages(c("sf","dplyr","readr","tidyr","matrixStats","MatchIt","prioritizr","highs"), repos="https://cloud.r-project.org")'
```

`sf` 需要系统库 GDAL/GEOS/PROJ；若安装失败，请用服务器上已有的 R 环境或 conda 环境。

---

## 3. 同步代码

在**本地**终端执行（把最新脚本推到服务器）：

```bash
rsync -avz --include='*.R' --exclude='*' \
  "$HOME/Documents/New project/bird_dynamic_occupancy_analysis/code_v3/" \
  dingchenchen@162.105.149.23:"~/Documents/New project/bird_dynamic_occupancy_analysis/code_v3/"
```

或直接在服务器上从 GitHub 拉取：

```bash
git clone https://github.com/dingchenchen6/china-bird-community-dynamics.git
cp china-bird-community-dynamics/code/code_v3/*.R ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3/
```

---

## 4. 一键运行

```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis
bash <本仓库>/server_run/run_mechanism_and_conservation.sh
```

驱动脚本会自动：检查 R 与依赖包 → 校验 psi 维度 → 校验保护区数据 → 按 30→31→32 顺序执行 → 每个脚本单独记日志 → 末尾汇总产出表。

**分步运行**（推荐首次这样做，便于定位问题）：

```bash
bash run_mechanism_and_conservation.sh 30     # 先跑机制归因
bash run_mechanism_and_conservation.sh 31     # 再跑保护成效
bash run_mechanism_and_conservation.sh 32     # 最后跑保护规划
```

**长时间任务建议放 tmux**：

```bash
tmux new -s analysis
bash run_mechanism_and_conservation.sh
# Ctrl+B 然后按 D 脱离；重连： tmux attach -t analysis
```

---

## 5. 三个脚本分别做什么、看什么结果

### 脚本 30：物种数增长的机制归因

| 输出表 | 关键判读 |
|---|---|
| `table_range_decomposition_summary_*.csv` | `pct_expansion` vs `pct_infilling` 谁占优：**扩张主导=气候位移特征；填充主导=种群恢复特征** |
| `table_colonization_extinction_*.csv` | 净增长来自"定殖多"还是"灭绝少" |
| `table_centroid_shift_*.csv` | 分布重心是否系统性**北移+上移**（气候变暖指纹） |
| `table_cti_*.csv` | 群落温度指数是否持续上升（热适应种占比上升） |
| `table_guild_contribution_*.csv` | 哪些生境类群贡献了物种数增量。**若高度集中于水鸟等少数类群，即解释了功能多样性为何没跟上** |

脚本还会自动校验一个恒等式：每格平均物种数 = (物种数/网格数) × 平均物种分布面积。屏幕会打印校验残差，应接近 0。

### 脚本 31：保护成效与保护空缺

| 输出表 | 关键判读 |
|---|---|
| `table_pa_grid_coverage_*.csv` | 各网格保护区覆盖率（按级别/类型/年代分层） |
| `table_pa_matching_att_*.csv` | **核心**：匹配后区内外趋势差异。若 `trend_fric_prob`/`trend_rao_q` 的 ATT 显著为正 = **保护区延缓了功能同质化** |
| `table_pa_type_effects_*.csv` | 湿地类 vs 森林类保护区效应对比（检验"湿地恢复驱动水鸟回归"假说） |
| `table_pa_did_*.csv` | 2000–2024 新建保护区的双重差分，正系数=保护后相对改善 |
| `table_conservation_gap_*.csv` | 各物种受保护比例 vs 30% 目标；若"栖息地宽度 vs 保护比例"正相关 = **特化种存在保护空缺** |

### 脚本 32：系统保护规划

| 输出表 | 关键判读 |
|---|---|
| `table_priority_solutions_*.csv` | 三套方案（物种数导向/功能独特性导向/同质化风险导向）逐网格选中结果 |
| `table_priority_overlap_*.csv` | **最重要**：方案间 Jaccard 重叠度。**重叠度低 = 只按物种数规划会系统性错过阻止同质化的关键区域** |
| `table_priority_pa_shortfall_*.csv` | 现有保护区对优先区的达成率；`n_priority_unprotected` = 扩建首选目标网格 |

---

## 6. 跑完之后

把这些表回传给用户（或提交到 GitHub 仓库的 `results/results_v3/`）：

```bash
# 在本地执行，拉回结果
rsync -avz --include='table_range_*' --include='table_colonization_*' \
  --include='table_centroid_*' --include='table_cti_*' --include='table_guild_*' \
  --include='table_pa_*' --include='table_conservation_gap_*' --include='table_priority_*' \
  --exclude='*' \
  dingchenchen@162.105.149.23:"~/Documents/New project/bird_dynamic_occupancy_analysis/results_v3/" \
  "$HOME/Desktop/analysis_results/"
```

---

## 7. 常见问题

| 症状 | 原因与处理 |
|---|---|
| `psi is 3D; need 4D` | psi 数组缺 period 维 → 重跑 `05_postprocess_diversity_extended.R` |
| `PA shapefile not found` | 保护区数据没上传 → 见 §2.2 |
| `MatchIt 未安装，跳过匹配` | 装包：`install.packages("MatchIt")` |
| `prioritizr 未安装` | 装包：`install.packages(c("prioritizr","highs"))` |
| `sf` 安装失败 | 缺 GDAL/GEOS/PROJ 系统库；改用 conda 环境或联系服务器管理员 |
| 脚本 31 报 `table_trend_summary_* not found` | 需先有 05 的扩展后处理产出 |
| 内存不足 | psi 数组较大；可减少 `PSI_MAX_DRAWS`（在 `00_config.R`）后重跑 05 |

---

## 8. 重要提醒（结果解释的边界）

1. **保护区矢量只覆盖 1028/3376 条**保护区（偏国家级与大面积），结论偏向大型保护区，写作时必须声明。
2. **匹配不能排除未观测混杂**，保护成效结论应表述为"与保护相关的差异"；**双重差分（脚本 31 的 DiD 部分）才是更强的准实验证据**。
3. 环境关联一律是**空间结构化的相关**，不是因果效应。
4. 脚本 30 的"扩张 vs 填充"分解依赖核心区阈值（默认 0.5），脚本已同时输出 0.3 与 0.7 的敏感性结果，报告时应说明结论是否随阈值改变。
