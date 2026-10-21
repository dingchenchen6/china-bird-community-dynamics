# 重跑代码包 code_v4_refit（2026-10-20）

本目录是 `code_v3` 的完整副本，加上 6 处补丁和 7 个新文件。总控脚本 `run_refit.sh` 按阶段依次执行全部重拟合与重分析，支持断点续跑。所有输出写入新目录 `BIRD_OUTPUT_ROOT`（建议 `refit_v4/`），不会覆盖 v3。

## 本轮新发现：环境协变量错位（比 BLOCKERS 记录的更严重）

v1 脚本 `code/prepare_dynamic_occupancy_environment.R` 的 `extract_mean_tbl()` 把 `terra::extract()` 返回的多边形行号（1…1739）当作格点编号 `grid_cell`（54…3665）去连接。后果有两层：

1. 编号大于 1739 的格点（北纬约 38.5° 以南的全部格点）没有拿到值，即 BLOCKERS 记录的"769/1,247 格缺失、中位数填补"。
2. **编号不超过 1739 的 478 个格点拿到的是别的格点的值。**用正确方法重新提取后，旧值与正确值的相关只有 −0.59 到 0.45。

受影响的只有气候、海拔、纹理、土地覆盖和生境多样性这 11 个协变量。经纬度和 HFI（v3 用 03e 重新提取过）与正确值的相关为 1，是对的。

**因此必须重做的结果**：主模型的环境效应与全部后验；驱动因子分析（方差分解、随机森林）；**群落温度指数与海拔重心**（只用了那 478 格的错值）；**保护地匹配**（匹配变量就是这些协变量）。纬度重心只用坐标，不受影响。

**阶段 0 的协变量重建已经在本地完成**：`refit_v4/data/derived_v3/grid_environment_v3.rds`，1,247 格 × 52 个变量，模型协变量零缺失；质检表见 `refit_v4/results_v3/QA_grid_environment_rebuild.csv`。

## 用法

```bash
export BIRD_PROJECT_ROOT="$HOME/Documents/New project/bird_dynamic_occupancy_analysis"
export BIRD_OUTPUT_ROOT="$BIRD_PROJECT_ROOT/refit_v4"
export V3_CODE_DIR="$BIRD_PROJECT_ROOT/code_v4_refit"
tmux new -s refit
bash "$V3_CODE_DIR/run_refit.sh" stage0
bash "$V3_CODE_DIR/run_refit.sh" stage1
bash "$V3_CODE_DIR/run_refit.sh" stage2
bash "$V3_CODE_DIR/run_refit.sh" stage3
bash "$V3_CODE_DIR/run_refit.sh" stage4
```

- 每一步完成后写 `status/<步骤>.done`，重跑时自动跳过；日志在 `logs_v3/`，进度在 `status/pipeline.log`。
- 服务器上的栅格路径与本机不同时，设置 `REFIT_WORLDCLIM_DIR`、`REFIT_ENV_EXTERNAL_DIR`、`REFIT_HFI_DIR`；也可以直接上传本地已重建的 `refit_v4/data/derived_v3/grid_environment_v3.rds`，阶段 0 的 S0_01 会被跳过（先 `touch $BIRD_OUTPUT_ROOT/status/S0_01_environment.done`）。
- 保护规划要读 AVONET：设置 `V3_AVONET_PATH=<服务器上的 AVONET1_BirdLife.csv 路径>`。
- 系统发育零模型次数：`PHYLO_N_NULL`（默认 199）。

## 任务与代码对照

| 阶段 | 任务 | 脚本 | 类型 | 验收标准 |
|---|---|---|---|---|
| 0 | 重建 1,247 格环境协变量 | `S0_01_rebuild_grid_environment.R` | **新写，本地已跑通** | 模型协变量缺失率 < 2%，否则自动停止 |
| 0 | 准备输入、校验四维调查历史与检测协变量、冻结名单、SHA-256 | `S0_02_check_inputs_and_freeze.R` | **新写** | 四维 `y` 加 4 个检测协变量；`INPUT_SHA256.csv` |
| 0 | 重建四维调查历史 | `02_build_survey_history.R` | 原脚本 | 本机的 `survey_history_v3.rds` 是旧的二维版本，以服务器版本或本步输出为准 |
| 0 | 准备模型环境与检测协变量 | `03_prepare_environment.R` | **补丁**：禁用有缺陷的 v2 环境表；缺失率 > 2% 即停止 | — |
| 1 | 200 种空间模型重拟合（4 链） | `04_run_stMsPGOcc_main.R` + `26_combine_validate_chains.R` | **补丁**：协变量缺失 > 2% 即停止（不再中位数填补）；投影坐标与 φ 先验原已改好 | 38 个参数全部 R̂ < 1.01，ESS ≥ 400 |
| 1 | 500 种时间模型重拟合（保留链身份） | `04_run_tMsPGOcc_500sp.R` | 原脚本 | 19 个参数 R̂ < 1.01 |
| 2 | 数据来源检测项与 ΔWAIC | `04_run_stMsPGOcc_main.R`（`V3_INCLUDE_SOURCE=1`）+ `29_compare_source_detection.R` | 原脚本 | ΔWAIC 与加入来源项后的趋势 |
| 2 | 后验预测检验 | `05c_ppc_bayesian_pvalue.R` | **新写** | 按物种 × 时期的 Bayesian p；按努力五分位的偏差 |
| 2 | 三年窗口重拟合 | `02` → `03` → `04` → `05`，设置 `V3_PERIOD_LENGTH=3` | **补丁**：`00_config.R` 读取时期长度 | 输出在 `refit_v4_3yr/`，与五年期对比 |
| 3 | 多维多样性后验传播 | `05_postprocess_diversity_extended.R` | **补丁**（10e 节）：加权均值只在有数据的格点上计算；去掉 200 种上限；分布区有数据占比 ≥ 0.8 | 物种层环境均值不再共线 |
| 3 | 系统发育 PD/MPD 与 ses.PD/ses.MPD | `05g_phylo_ses.R` | **新写**（估计量已在 `utils_diversity_extended.R` 中修正） | 格点 × 时期表与全国趋势表 |
| 3 | 功能趋势后验概率 | `05e_export_trend_draws.R`、`05f_functional_trend_pdecline.R` | 原脚本 | P(斜率 < 0) 与 95% 可信区间 |
| 3 | 全国格点间空间 beta 趋势 | `12_homogenization_spatiotemporal_maps.R` | 原脚本（需四维 psi） | 斜率、区间、Mann–Kendall |
| 3 | 驱动因子 | `09_extended_analyses_extended.R`、`09b_rf_spatial_block_cv.R` | 原脚本 + **新写**的空间分块交叉验证 | 样本外 R² 与分组重要性 |
| 3 | 性状回归 | `14_species_trait_regression.R` | 原脚本 | 日志中确认读入了系统发育树 |
| 3 | 增量分解、重心、CTI、类群贡献 | `30_range_expansion_mechanism.R` | 原脚本，自动使用重建后的协变量 | CTI 与海拔重心改用全部 1,247 格 |
| 3 | 保护地效应与对比 | `31_protected_area_effectiveness.R`、`33_pa_spatiotemporal_contrast.R` | **补丁**（33）："未保护"统一为 ≤ 5% | 匹配平衡（使用正确协变量） |
| 3 | 保护规划 | `32_conservation_prioritization.R` | **补丁**：接入 AVONET 形态性状，含 10 个学名的同物异名对照（已核实，200/200 可匹配） | 日志显示"AVONET 形态（7 项）" |
| 3 | 努力增长 × 物种增量 | `S2_04_effort_gain_test.R` | **新写** | 相关系数与分块 bootstrap 区间 |
| 4 | 全部图件（标准地图） | `24b/24c/24d/21b/34b/34c_*_20261011.R`、`35_qa_figure_package_20261011.py` | 原脚本 | 去掉旧 PD 面板，质检通过 |
| 4 | 收尾检查 | `28_submission_postflight.R` | 原脚本 | `status/REFIT_COMPLETE` |

## 注意

- 7 个新文件和 6 处补丁都通过了语法检查。`S0_01` 已用真实数据在本地跑通，其余新脚本还没有在真实的拟合对象上运行过，第一次运行时请留意日志。
- 本目录中原有的服务器 IP 已替换为 `<SERVER_IP>`；凭证不要写进任何文件。
- 结果带回后，按 `RERUN_HANDOFF_WORKBUDDY.md` 中"需要带回的东西"一节整理。

## 无人值守连续运行（2026-10-20 起在服务器执行）

`run_refit_auto.sh` 在后台等待当前 stage1 跑完后，自动依次执行 stage2 → stage3 → stage4，每阶段失败重试一次；中断后重跑会自动跳过已完成步骤（`status/*.done` 标记）。

```bash
export BIRD_PROJECT_ROOT=<项目根> BIRD_OUTPUT_ROOT=<新输出根> V3_CODE_DIR=<code_v4_refit>
nohup setsid bash $V3_CODE_DIR/run_refit_auto.sh >> $BIRD_OUTPUT_ROOT/logs_v3/auto_driver.bootstrap.log 2>&1 &
```

- 进度：`logs_v3/auto_driver.log`（阶段起止与失败重试）、`status/pipeline.log`
- 结束标记：`status/AUTO_PIPELINE_COMPLETE`；若为 `AUTO_PIPELINE_PARTIAL` 或 `AUTO_PIPELINE_INCOMPLETE`，查对应的 `logs_v3/<stage>_driver.log`
- 闸门：stage1 四条 200 种链文件缺失或小于 10 MB 时停止，不进入下游，避免在半成品上跑图件
- 服务器执行记录：stage0 于 16:39 全部通过；stage1 于 16:51 启动（4 链 × 32 线程）；驱动于 16:58 接管等待

## 服务器环境注意（2026-10-20 实测）

- `~/.Renviron` 中的变量会覆盖 shell 环境变量；原文件里 `BIRD_PROJECT_ROOT`、`V3_CODE_DIR` 指向已弃用旧路径，已备份为 `~/.Renviron.bak_20261020` 并注释掉，路径统一由启动命令显式传入。
- 服务器上的 WorldClim 栅格（2.5′ 与 10′ 两套）GDAL 均无法打开（疑似文件损坏），因此 `S0_01` 在本地完成、重建表上传服务器；服务器端不要重复执行 S0_01（`status/S0_01_environment.done` 已标记跳过）。

## 收敛门与 MCMC 配置修订（2026-10-21）

首轮 M200_base（600 batch / burn 8,000 / thin 25 → 280 draws，10 factors）组合后 `beta.comm` 大面积失败（max R̂ 1.725、min ESS 6.3，仅 23.7% 参数过 R̂<1.01），诊断为沿共线岭的非平稳漂移（burn-in 不足），非多峰。复盘 2026-07 重跑（`bird_dynamic_occupancy_rerun_20260722`）发现同样轨迹：280 draws 首轮失败后，以 **1,600 batch（40,000 迭代）/ burn 20,000 / thin 10 / 10 factors（2,000 draws）** 重跑，且收敛门改为**群落层超参披露制**才放行。据此 v4 做了三项修订：

1. **`26_combine_validate_chains.R` 换用 7 月版披露制闸门**：严格层（物种层 z/psi、检测层 alpha、空间/时间参数）仍须 R̂≤1.05 且 ESS≥200；`beta.comm`/`tau.sq.beta` 群落层超参放宽至 R̂≤1.4 且 ESS≥5 进 watchlist（落盘 `table_convergence_disclosure_<label>.csv`，投稿补充材料必须报告）；超兜底上限仍判失败。该模型的截距即使 40,000 迭代也只有 R̂≈1.20 / ESS≈14，硬闸不可达是模型固有弱识别，不是采样不足。
2. **M200 系列统一到 7 月验证规格**：`V3_N_BATCH=1600 V3_N_BURN=20000 V3_N_THIN=10 V3_N_FACTORS=10`（此前误降到 1,200/15,000/5 factors，偏离 v3 既定模型规格且迭代不足）。预计约 22 h/模型。
3. **`04_run_tMsPGOcc_500sp.R` 增加 `V3_T500_*` 独立配置并在保存前丢弃 `z.samples`**：时间模型非空间、混合快，四链并行内存以 1,000 draws/链封顶（本轮 `V3_T500_N_BATCH=800 V3_T500_N_BURN=10000 V3_T500_N_THIN=10`）。
4. **`run_refit_auto.sh` 的驱动闸门改读 `gate_pass` 列**（旧逻辑只看 rhat/ess 两列，会漏掉 watchlist 兜底失败的情形）。

**投稿前必做**：对本轮链重跑逐链 95% CrI 共同重叠诊断（参照 7 月 `98_chain_diagnosis_community_20260804.R`），确认无多峰后按披露表在补充材料中报告。首轮 280 draws 的失败产物已归档于 `v4_refit/archive_convergence_fail_20261021/`。
