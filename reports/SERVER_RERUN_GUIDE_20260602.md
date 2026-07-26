# Server 补算指南 —— 闭环 GCB 稿件 3 个 [PENDING]（2026-06-02）

> 在 **23 服务器**（<SERVER_IP>，有 v3 的 fit / psi / metric checkpoint）上运行。
> 目标：把 `MANUSCRIPT_GCB_v6` 里剩余 3 个 `[PENDING]` 用实测值闭环。
> 三项相互独立，可并行；难度：项 2 最易，项 1 中，项 3 最重。

## 通用环境变量
```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3
export V3_CODE_DIR="$PWD"
export V3_RUN_LABEL="v3_full_200sp_ar1_spatial"   # 与现有产物一致
export V3_OMP_THREADS=8                            # 按机器核数
```

---

## 项 1 — 空间 beta（同质化金标准证据）｜脚本已修复
`12_homogenization_spatiotemporal_maps.R` 第 92–97 行已修：3D psi 不再静默产 NA，改为显式报错；4D psi 正常计算。**只需重跑**（用当前 4D `psi_samples_thinned_<run_label>.rds`）：
```bash
Rscript 12_homogenization_spatiotemporal_maps.R
```
- 预期输出：`results_v3/table_homogenization_trend_v3_full_200sp_ar1_spatial.csv` 的 `mean_sorensen` 不再是 NA，5 期给出跨格点相似度 + Mann-Kendall 趋势。
- 若仍报 "psi is 3D..." 错误：说明该 run 的 `psi_samples_thinned` 仍是旧 3D 版本 → 先重跑 `05_postprocess_diversity_extended.R`（其 FIX #3 会保存 4D psi + n_periods），再跑 12。
- **回填稿件**：3.3 节 `[PENDING: spatial β slope ...]` → 填 Sørensen 各期均值下降的斜率/MK 检验 P 值。

## 项 2 — 功能指标 P(decline)（闭环 homogenization 主张）｜可直接跑
```bash
Rscript 05f_functional_trend_pdecline.R
```
- 复用 `metric_arrays_checkpoint_<run_label>_extended.rds` 与 `theil_sen_slope`，对每个后验 draw 算"全国均值时间序列"的斜率 → 后验分布。
- 预期输出：`results_v3/table_functional_trend_pdecline_<run_label>.csv`，列含 `slope_mean / slope_l95 / slope_u95 / P_decline`；控制台直接打印 trait_volume、rao_q 的可粘贴措辞。
- **回填稿件**：3.2 节 `[trait-volume slope ... 95% CrI PENDING, P(decline) PENDING; Rao's Q ...]` → 填实测 CrI 与 P(decline)。
  - 若 `P_decline ≥ 0.95`（CrI 不含 0）→ 用强措辞 "functional homogenization"；
  - 若 CrI 跨 0 → 降级为 "taxonomic expansion without detectable functional expansion"（稿件已内置这个开关）。

## 项 3 — source-term ΔWAIC（混淆控制最后一环）｜最重，需先对齐 2 处
`04d_source_detection_refit.R` 会读 base fit、构造来源协变量、重跑并比 WAIC。**运行前务必核对脚本内两处 `## >>> ADJUST`**：
- **ADJUST #1**：`survey$records`（replicate 级长表）需含 `site_index / block_id / year_in_block / source`。若你的 survey history 未在 replicate 级保留 `source`，需先在 `02_build_survey_history.R` 增加该列（source 在 01_merge 的 events 里已有：birdwatch/ebird/gbif）。
- **ADJUST #2**：把脚本里的 `inits / priors / tuning / n_batch` 替换为 `04_run_stMsPGOcc_main.R` 第 357–405 行的**完全相同**设定（否则 WAIC 不可比）。
```bash
Rscript 04d_source_detection_refit.R
```
- 预期输出：`results_v3/table_source_detection_waic_<run_label>.csv`（base vs +source 的 WAIC、ΔWAIC）+ 新 fit `stMsPGOcc_fit_<run_label>_srcdet.rds`。
- **回填稿件**：3.4 节 `[PENDING: ΔWAIC and source-term trend slope]` → 填 ΔWAIC；并核对加入 source 项后群落趋势斜率是否基本不变（稳健性）。

---

## 完成后
把三张结果表（`table_homogenization_trend_*`、`table_functional_trend_pdecline_*`、`table_source_detection_waic_*`）发我，我据实测值：
1. 精确回填 GCB 稿件全部 `[PENDING]`；
2. 依 P(decline) 锁定 "functional homogenization" 措辞档位；
3. 重新生成 DOCX。
