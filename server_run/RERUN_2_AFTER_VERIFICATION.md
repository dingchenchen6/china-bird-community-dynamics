# 第二轮重跑说明（核查修复后）

> 面向：Kimi 或其他助手，在计算服务器上执行。**本文自包含**，无需先前对话上下文。
> 背景：首轮产出中有三个结论经核查不成立（详见 `reports/RESULTS_VERIFICATION_20260728.md`）。
> 脚本已修复，需重跑以取代旧结果。

---

## 0. 一句话说明要做什么

重跑脚本 **30 / 32 / 33 / 34**（31 无需重跑），用修复后的代码取代首轮的错误结果。

---

## 1. 首轮哪些结论被推翻（重跑前务必知情）

| 首轮结论 | 核查结果 | 处理 |
|---|---|---|
| 「保护年限越久功能维度越好」p=0.0007 | ❌ 伪重复导致的假显著。直接重算 r = **+0.082**，无实质关系 | 模型加随机斜率；新增 `effect_size_r` |
| 「period×保护交互 t = −14.6」 | ❌ 同一问题，p 值不可信；**但效应方向由原始均值独立支持** | 模型加随机斜率；改报效应量 |
| 「物种数导向 ≈ 功能导向规划 J=0.96」 | ❌ 构造性假象：两方案只差物种层权重，其 CV 仅 0.15 | 方案 B 重建为站点层功能目标 |
| 「保护地内同质化更慢」 | ⚠️ 方向成立但 `sd=0` 无法检验 | 改 bootstrap + 置换检验 |
| 「气候驱动被排除」 | ✅ 成立（重心无位移、CTI 下降） | 保留 |
| 「湿地类群主导增量」 | ✅ 成立（21.5% 物种贡献 38.0% 增量） | 保留，需脚本 30 重跑固化 |

**重要**：重跑后若 `effect_size_r` 仍 |r| < 0.1，说明剂量-反应确实不存在——**这是正确结果，不是失败**，不要试图调参使其"显著"。

---

## 2. 环境与前置（与第一轮相同）

```bash
ssh dingchenchen@162.105.149.23      # 密码由用户单独提供
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis
export V3_CODE_DIR="$PWD/code_v3"
export V3_RUN_LABEL="v3_full_200sp_ar1_spatial"
export V3_OMP_THREADS=8
```

前置条件（应已满足）：4 维 `psi_samples_thinned`、保护区数据、R 包。本轮新增依赖：**`lme4`**（随机斜率模型必需）。

```bash
Rscript -e 'if (!requireNamespace("lme4", quietly=TRUE)) install.packages("lme4", repos="https://cloud.r-project.org")'
```

---

## 3. 同步修复后的代码

```bash
cd ~ && rm -rf china-bird-community-dynamics
git clone https://github.com/dingchenchen6/china-bird-community-dynamics.git
cp china-bird-community-dynamics/code/code_v3/*.R \
   ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3/
```

确认拿到的是修复版（应能搜到随机斜率写法）：

```bash
grep -c "period_num | grid_cell" ~/Documents/New\ project/bird_dynamic_occupancy_analysis/code_v3/33_pa_spatiotemporal_contrast.R
# 期望输出 >= 2；若为 0 说明代码没更新
```

---

## 4. 备份首轮结果，再重跑

```bash
cd ~/Documents/New\ project/bird_dynamic_occupancy_analysis
mkdir -p results_v3_round1_backup
cp results_v3/table_pa_*.csv* results_v3/table_priority_*.csv* \
   results_v3/table_range_*.csv* results_v3/table_guild_*.csv* \
   results_v3_round1_backup/ 2>/dev/null

bash ~/china-bird-community-dynamics/server_run/run_mechanism_and_conservation.sh 30 32 33 34
```

**不必重跑 31**（匹配与 DiD 逻辑未变；其产出 `table_pa_grid_coverage` 是 33 的输入）。
若 31 的覆盖率表不存在，则需先跑 31。

建议放 tmux：`tmux new -s rerun2`，Ctrl+B 再 D 脱离。

---

## 5. 本轮新增/变更的产出

| 表 | 变化 | 判读 |
|---|---|---|
| `table_pa_age_dose_response` | **新增 `effect_size_r` 列** | **以 r 为准**：\|r\| < 0.1 = 无实质剂量-反应，不论 p 多小 |
| `table_pa_interaction_tests` | 新增 `n_grid`、`ranef` 列 | `ranef` 应为 `random_slope`；若显示 `intercept_only` 说明模型未收敛，p 值仍不可信 |
| `table_pa_beta_contrast` | **列名全变**：`*_lwr/_upr` 取代 `*_sd`，新增 `beta_gap_obs`、`perm_p` | 仅当 `beta_gap` 区间不含 0 **且** `perm_p < 0.05`，才可称保护区抑制趋同 |
| `table_pa_beta_rate` | **全新** | difference 行的区间不含 0 才可断言内外同质化速率不同 |
| `table_priority_overlap` | 方案 B 已重建 | 与首轮 J=0.96 不可比；重跑值才有效 |
| `table_guild_contribution` | 首轮未生成（路径 bug 已修） | 应显示湿地类群超额贡献 |
| `fig_pa5b_homogenization_rate_*` | **全新图** | 速率差异的区间图 |

运行时留意屏幕上这几行提示：

```
[32] 物种层功能权重 w_func 的 CV = ...     ← 若仍 < 0.2，方案 B 仍主要靠站点层区分
[32] 功能独特性性状空间：AVONET 形态 / 生态属性评分
[33] 注：有效样本量为网格数 n_grid，而非观测数 n_obs
[33] >>> 判读以 effect_size_r 为准
```

---

## 6. 跑完回传

```bash
# 本地执行
rsync -avz --include='table_pa_*' --include='table_priority_*' \
  --include='table_guild_*' --include='table_range_*' \
  --include='fig_pa*' --exclude='*' \
  dingchenchen@<SERVER_IP>:"~/Documents/New project/bird_dynamic_occupancy_analysis/results_v3/" \
  ~/Desktop/rerun2_results/
```

⚠️ 服务器上部分表名带 `.csv.csv` 双后缀（31/32 脚本写法所致），回传后需改名：

```bash
cd ~/Desktop/rerun2_results && for f in *.csv.csv; do mv "$f" "${f%.csv}"; done
```

---

## 7. 判读原则（防止再次误读）

1. **效应量优先于 p 值**。本项目所有面板/多期数据的 p 值都可能被伪重复放大；有 `effect_size_r` 就以它为准。
2. **模型系数与原始均值冲突时，信原始均值**——首轮的错误正是这样暴露的。
3. **阴性结果是结果**。剂量-反应不存在、保护成效不显著，都是可发表的诚实结论，不要调参凑显著。
4. **注意天花板效应**：高保护网格功能体积基线本就更高（1.4105 vs 1.3794），降幅更大可能是回归均值，不等于"保护有害"。
5. 保护区边界仅覆盖 1028/3376 条且偏大型国家级，**成效估计偏保守**。
