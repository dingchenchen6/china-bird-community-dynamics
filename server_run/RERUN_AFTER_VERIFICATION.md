# 核查后的重跑指令（供 Kimi 在服务器执行）

> 首轮结果经复核发现三处问题（详见 `reports/RESULTS_VERIFICATION_20260728.md`），
> 脚本 30/32/33/34 已修复。本文件是**自包含**的重跑说明。
> 目标：用修正后的方法重算，并生成两张新表 + 一张新图。

---

## 0. 为什么要重跑（一分钟背景）

| 首轮结论 | 问题 | 修复 |
|---|---|---|
| "保护年限越久功能越好"（p=0.0007） | 混合模型只有随机截距，保护年限是网格层常量 → 每格 5 期被当成独立信息，有效样本 489 虚增到 2445，SE 被低估 | 改随机斜率 `(period_num \| grid_cell)`，并输出不依赖模型的 `effect_size_r` |
| "period×保护交互 t=−14.6" | 同一问题，n_obs 6235 实为 1247 格×5 期 | 同上 |
| "物种数导向≈功能导向规划 J=0.96" | 两方案只差物种层权重，其 CV 仅 0.15，站点效用加总后排序几乎不变；而方案 C 用站点层权重，结构不对等 | 方案 B 增加站点层功能目标；性状空间优先 AVONET 形态 |
| 空间 beta `sd = 0`，无法检验 | 保护组恰 91 格、无放回抽样每次抽全 → 方差恒为 0 | 改**有放回 bootstrap + 置换检验**，并新增同质化**速率**对比 |
| `table_guild_contribution` 未生成 | 脚本 30 把 `run_label` 误传为 `v3_file` 的扩展名参数 | 路径已修正 |

---

## 1. 前置检查

```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis
export V3_CODE_DIR="$PWD/code_v3"
export V3_RUN_LABEL="v3_full_200sp_ar1_spatial"

# psi 必须是 4 维
Rscript -e 'x <- readRDS("data/derived_v3/psi_samples_thinned_v3_full_200sp_ar1_spatial.rds"); print(dim(x$psi_samples_thinned))'
# R 包
Rscript -e 'for (p in c("lme4","prioritizr","highs","sf","dplyr","readr","tidyr","matrixStats","ggplot2","patchwork","scico")) if (!requireNamespace(p, quietly=TRUE)) cat("缺:", p, "\n")'
```

## 2. 同步修复后的代码

```bash
cd ~ && git clone https://github.com/dingchenchen6/china-bird-community-dynamics.git repo-latest 2>/dev/null || (cd repo-latest && git pull)
cp ~/repo-latest/code/code_v3/{30,32,33,34}_*.R \
   ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3/
```

## 3. 重跑（按顺序，31 无需重跑）

```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis
bash ~/repo-latest/server_run/run_mechanism_and_conservation.sh 30
bash ~/repo-latest/server_run/run_mechanism_and_conservation.sh 32
bash ~/repo-latest/server_run/run_mechanism_and_conservation.sh 33
bash ~/repo-latest/server_run/run_mechanism_and_conservation.sh 34
```

脚本 33 的 bootstrap（B=200）+ 置换检验（500 次）× 5 期较慢，建议放 tmux：

```bash
tmux new -s rerun
# 跑完 Ctrl+B 再按 D 脱离；重连 tmux attach -t rerun
```

---

## 4. 重点核对这些输出

### 新增表 1：`table_pa_beta_rate_*.csv`
同质化**速率**的内外对比（bootstrap 95% 区间，B=200）：

| 判读 | 含义 |
|---|---|
| `difference` 行的 `lwr`–`upr` **不含 0** | 保护地内外同质化速率确有差异，可作正式结论 |
| **含 0** | 差异不能与偶然区分，只能作描述性证据 |

### 新增表 2：`table_guild_contribution_*.csv`
生境类群对物种数增量的贡献。**预期**（本地已用同一数据先行验算）：湿地类群占物种数 21.5% 却贡献 **38.0%** 的增量（1.77 倍超额），森林 29% 物种仅贡献 18.6%。**若服务器结果与此明显不符，请回报——说明存在数据版本差异。**

### 已有表的关键变化

| 表 | 看什么 |
|---|---|
| `table_pa_age_dose_response` | **以 `effect_size_r` 为准，不看 p 值**。\|r\| < 0.1 即判定无实质剂量-反应（首轮直接重算为 r = 0.08） |
| `table_pa_interaction_tests` | 新增 `n_grid`（有效样本）与 `ranef` 列。若 `ranef = intercept_only` 说明随机斜率未收敛，该行结果需谨慎 |
| `table_pa_beta_contrast` | 新增 `beta_gap_lwr/upr` 与 `perm_p`；仅当区间不含 0 **且** `perm_p < 0.05` 才可称保护抑制了趋同 |
| `table_priority_overlap` | 方案 B 已重建。运行日志会打印 `w_func` 与站点层的 CV——**若站点层 CV 仍很小，说明功能目标依然无力，需回报** |

### 新增图
`fig_pa5b_homogenization_rate_*`：同质化速率差异的 bootstrap 区间图。

---

## 5. 回传结果

```bash
# 本地执行
rsync -avz --include='table_pa_beta_rate*' --include='table_guild_contribution*' \
  --include='table_pa_age_dose_response*' --include='table_pa_interaction_tests*' \
  --include='table_pa_beta_contrast*' --include='table_priority_*' \
  --include='fig_pa5b*' --exclude='*' \
  dingchenchen@<SERVER_IP>:"~/Documents/New project/bird_dynamic_occupancy_analysis/results_v3/" \
  "$HOME/Desktop/rerun_results/"
```

---

## 6. 判读原则（避免重蹈覆辙）

1. **效应量优先于 p 值**。本轮出问题的两个结论都是 p 值极小但效应量近零。
2. **有效样本量是网格数，不是观测数**。1247 格 × 5 期 ≠ 6235 个独立观测。
3. **回归系数与原始分组均值方向不一致时，以原始均值为准**，并检查模型设定。
4. **方案对比要结构对等**。用物种层权重 vs 站点层权重比较，差异可能来自构造而非生态。
5. **零方差是危险信号**。若某组的 bootstrap 标准差为 0，多半是抽样设计退化，不是真的没有不确定性。
