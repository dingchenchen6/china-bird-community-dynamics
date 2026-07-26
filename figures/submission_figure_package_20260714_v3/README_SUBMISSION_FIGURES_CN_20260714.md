# 中国鸟类群落动态研究：顶刊投稿图件包 v3

本目录是面向 Nature 类及其他高影响力期刊的正式主图交付包。图件按 183 mm 双栏宽度设计，最大高度不超过 170 mm；正文和刻度使用 Helvetica 5-7 pt，分面标签使用 8 pt 粗体小写字母。配色避免以红绿作为唯一编码，统计图保留轴线与刻度并取消背景网格。

## 核心内容

- `figures/`：6 幅投稿主图。
- `source_data/`：每幅图对应的 Source Data CSV。
- `qa/FIGURE_QA_CONTRACT_20260714.csv`：逐图结论、证据层级、尺寸、风险与导出契约。
- `FIGURE_LEGENDS_BILINGUAL_20260714.md`：可直接进入稿件的中英文自足式图注。
- `NATURE_SUBMISSION_MAIN_FIGURES_ATLAS_20260714.pptx`：6 幅主图的稳定浏览版。
- `SUBMISSION_FIGURE_MANIFEST_20260714.csv`：全部格式、尺寸、可编辑性和 Source Data 路径。

## 六幅主图

1. **Fig. 1 | 推断架构与证据边界**：区分 200 物种空间主模型和 500 物种时间广度扩展能够支持的结论。
2. **Fig. 2 | 五时期空间多维多样性**：200 物种空间模型的 6 个核心指标、5 个时期全国格局。
3. **Fig. 3 | 多样性变化的空间指纹**：分类、功能和系统发育趋势的空间解耦。
4. **Fig. 4 | 多维度变化解耦**：跨模型端点效应与 500 物种核心指标轨迹。
5. **Fig. 5 | 周转崩塌、赢家与性状**：晚期 nestedness、扩张主导和广栖息地物种获益。
6. **Fig. 6 | 检测偏差与空间驱动背景**：朴素趋势偏差、方向翻转、随机森林和变异分解。

## 导出格式

- 所有主图：PDF、PNG、450 dpi LZW TIFF、PPTX。
- 非地图图：另有 SVG；独立 PPTX 为原生 DrawingML，可编辑文字、点、线、坐标轴和图例。
- 地图图：PDF 保留矢量几何；PPTX 使用高分辨率 PNG 嵌入，避免数万多边形导致 PowerPoint 不稳定。
- 合集 PPTX：用于浏览和汇报，采用整图嵌入保证版式不漂移。

## 地图与推断边界

- 所有地图使用 `data/中国shp/省.shp`，南海区域直接纳入完整主图，不使用鹰眼图。
- Fig. 2-3 来自 200 物种空间主模型，可支持空间推断。
- 500 物种模型为非空间 `tMsPGOcc`，只承担时间广度与物种池一致性检验。
- 性状分析有效样本为完整案例 n=200。
- 驱动分析是观察性关联或预测重要性，不是因果效应。

## 全量增强版补充图谱

- `figures_top_journal_20260714_v2/`：11 幅跨模型综合图，均有 PNG、PDF、PPTX及图册。
- `figures_500sp_all_analysis_20260714_v2/`：24 幅 500 物种完整分析图，均有 PNG、PDF、PPTX及图册。
- 全量图谱已同步采用无网格、色盲友好语义色和统一 Helvetica 视觉系统。
- 旧版 `pd_prob` 和 `mpd_prob` 为全缺失值，在全量图谱中保持空白，不被误画为零变化。

## 复现命令

```bash
Rscript code_v3/21_make_top_journal_figures.R
Rscript code_v3/22_make_500sp_all_analysis_figures.R
Rscript code_v3/24_make_nature_submission_figures.R
```

正式投稿建议优先提交各图独立 PDF；PNG 和 TIFF 用于预览或期刊在线系统，PPTX用于作者内部编辑。
