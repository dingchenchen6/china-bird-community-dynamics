#!/usr/bin/env Rscript
## 20_concept_workflow_diagrams.R
##
## 1) 概念图：多样性三维（Taxonomic / Phylogenetic / Functional）+ 驱动 -> 占域 -> 群落动态
##    用 ggforce::geom_circle 画 Venn-like 三圆 + ggplot 注释
## 2) 技术流程图：raw 数据 → tMsPGOcc → psi 后验 → 指标计算 → 驱动归因 → 出版
##    用 DiagrammeR + DiagrammeRsvg + rsvg → PDF + PNG

suppressPackageStartupMessages({
  library(ggplot2); library(ggforce); library(dplyr); library(tibble); library(scales)
  library(DiagrammeR); library(DiagrammeRsvg); library(rsvg)
})
CODE_V2 <- Sys.getenv("V2_CODE_DIR",
  "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/code_v2")
source(file.path(CODE_V2, "utils_paths.R"))
source(file.path(CODE_V2, "utils_mapping.R"))
P <- ensure_v2_dirs()

message("[stage-20] concept + workflow diagrams")

## --- 1. 概念图 ----------------------------------------------------------

circles <- tibble(
  dim  = c("Taxonomic", "Phylogenetic", "Functional"),
  x    = c(-1.3,  1.3, 0),
  y    = c(-0.7, -0.7, 1.0),
  r    = c(2.2, 2.2, 2.2),
  fill = c("#0E5A78", "#6A4C93", "#C9784A"))

indices <- tibble(
  dim   = c(rep("Taxonomic", 3), rep("Phylogenetic", 5), rep("Functional", 6)),
  label = c("Richness","Shannon","Inverse Simpson",
            "Faith's PD","MPD","MNTD","NRI","NTI",
            "FRic","FEve","FDiv","FDis","Rao Q","CWM/CWSD"),
  x = c(-3.0,-3.0,-3.0,  2.4,2.4,2.4,2.4,2.4,  -0.5, 0.4, -0.5, 0.4, -0.5, 0.4),
  y = c( 0.4, 0.0,-0.4,  0.4,0.1,-0.2,-0.5,-0.8, 2.6, 2.35, 2.1, 1.85, 1.6, 1.35))

p_concept <- ggplot() +
  ggforce::geom_circle(data = circles,
                        aes(x0 = x, y0 = y, r = r, fill = fill),
                        alpha = 0.18, colour = NA) +
  ggforce::geom_circle(data = circles,
                        aes(x0 = x, y0 = y, r = r, colour = fill),
                        fill = NA, linewidth = 0.6) +
  geom_text(data = circles,
            aes(x = x, y = y + r - 0.18, label = dim, colour = fill),
            size = 5.2, fontface = "bold") +
  geom_text(data = indices,
            aes(x = x, y = y, label = label,
                colour = circles$fill[match(dim, circles$dim)]),
            size = 3.3, hjust = 0) +
  scale_fill_identity() + scale_colour_identity() +
  annotate("text", x = 0, y = -2.3,
           label = "Multidimensional community structure",
           size = 6, fontface = "bold", colour = "#1F1F1F") +
  annotate("text", x = 0, y = -2.65,
           label = "Each axis captures complementary aspects of biodiversity",
           size = 3.6, colour = "grey25") +
  coord_fixed(xlim = c(-4.2, 4.2), ylim = c(-3.2, 4.6)) +
  theme_void(base_family = "Helvetica") +
  theme(plot.margin = margin(8, 8, 8, 8),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                    margin = margin(b = 6)),
        plot.subtitle = element_text(size = 11, hjust = 0.5,
                                        colour = "grey25",
                                        margin = margin(b = 8))) +
  labs(title = "Multidimensional diversity of bird communities",
       subtitle = "Three complementary axes share information yet capture distinct ecological signal")

# 驱动 -> psi 链顶部
driver_y <- 3.7
drv_box <- tibble(
  xmin = c(-3.7, -1.7, 0.3, 2.3),
  xmax = c(-1.9, -0.3, 1.7, 3.7),
  ymin = driver_y - 0.4, ymax = driver_y + 0.4,
  label = c("Climate", "Land use", "Topography", "Human pressure"))
p_concept <- p_concept +
  geom_rect(data = drv_box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#F2E8D8", colour = "#7C7C7C", linewidth = 0.45) +
  geom_text(data = drv_box,
            aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = label),
            size = 3.6, colour = "#1F1F1F", fontface = "bold") +
  annotate("segment", x = 0, y = driver_y - 0.45, xend = 0, yend = 2.05,
           colour = "grey40", linewidth = 0.45,
           arrow = grid::arrow(length = grid::unit(0.20, "cm"),
                                type = "closed")) +
  annotate("text", x = 0.25, y = 2.85,
           label = "drives via\noccupancy psi", size = 3.3, colour = "grey25",
           hjust = 0)

save_dual_v3(p_concept, "fig_v3_concept_diversity", width = 11, height = 9)

## --- 2. 工作流 DAG -----------------------------------------------------

dot <- '
digraph workflow {
  graph [layout = dot, rankdir = LR, fontname = Helvetica, bgcolor = white,
         pad = 0.4, nodesep = 0.4, ranksep = 0.55]
  node  [fontname = Helvetica, shape = box, style = "rounded,filled",
         fontsize = 11, penwidth = 0.6, color = "#333333"]
  edge  [fontname = Helvetica, fontsize = 9, color = "#666666", penwidth = 0.8,
         arrowsize = 0.7]

  subgraph cluster_data {
    label = "Data integration"; color = "#0E5A78"; style = "rounded,dashed"
    eBird [label = "eBird / GBIF\nchecklist", fillcolor = "#E6F0F4"]
    BW    [label = "China Birdwatch\nPlatform",   fillcolor = "#E6F0F4"]
    iNat  [label = "iNaturalist",  fillcolor = "#E6F0F4"]
    Dedup [label = "Deduplication\nspecies x date x location\nx observer", fillcolor = "#F0DCC0"]
  }
  subgraph cluster_grid {
    label = "Spatial / temporal frame"; color = "#6A4C93"; style = "rounded,dashed"
    Grid    [label = "100 km grid\n(1308 cells)", fillcolor = "#E9DEF2"]
    Blocks  [label = "5-year primary\nperiods 2000-2024", fillcolor = "#E9DEF2"]
    Cov     [label = "Detection covariates\nevents / observers / duration", fillcolor = "#E9DEF2"]
  }
  subgraph cluster_model {
    label = "Bayesian model"; color = "#C9784A"; style = "rounded,dashed"
    Fit  [label = "spOccupancy::tMsPGOcc\n200 sp, 4 chains, AR1",
          fillcolor = "#F8E0D2", fontname = "Helvetica-Bold"]
    Diag [label = "MCMC diagnostics\n(R-hat, ESS, PPC)", fillcolor = "#F8E0D2"]
    Psi  [label = "Posterior psi\n200 thinned draws", fillcolor = "#F8E0D2"]
  }
  subgraph cluster_div {
    label = "Multidimensional diversity"; color = "#3C8C5A"; style = "rounded,dashed"
    Tax  [label = "Taxonomic\nrichness, Shannon", fillcolor = "#DCEFE2"]
    Phy  [label = "Phylogenetic\nPD, MPD, NRI, NTI\nMcTavish tree (clootl)", fillcolor = "#DCEFE2"]
    Fun  [label = "Functional\nFRic, FEve, FDiv\nFDis, Rao Q\nfundiversity", fillcolor = "#DCEFE2"]
    Beta [label = "Temporal beta\nBaselga + homogenisation", fillcolor = "#DCEFE2"]
  }
  subgraph cluster_inf {
    label = "Driver inference"; color = "#8B2E1E"; style = "rounded,dashed"
    VIF  [label = "VIF + correlation\nscreening", fillcolor = "#F5E0E0"]
    Brms [label = "brms + cmdstanr\ngp(lon,lat) + interactions", fillcolor = "#F5E0E0"]
    Imp  [label = "Variance decomp.\nLMG + RF + varpart", fillcolor = "#F5E0E0"]
  }
  subgraph cluster_out {
    label = "Outputs"; color = "#1F1F1F"; style = "rounded,dashed"
    Maps   [label = "Spatial maps\ntime slices + trends", fillcolor = "#EFEFEF"]
    Curves [label = "Time series + CRI", fillcolor = "#EFEFEF"]
    Paper  [label = "Manuscript +\nPPT decks", fillcolor = "#EFEFEF"]
  }

  eBird -> Dedup ; BW -> Dedup ; iNat -> Dedup
  Dedup -> Grid ; Dedup -> Blocks ; Dedup -> Cov
  Grid -> Fit ; Blocks -> Fit ; Cov -> Fit
  Fit -> Diag ; Fit -> Psi
  Psi -> Tax ; Psi -> Phy ; Psi -> Fun ; Psi -> Beta
  Tax -> VIF ; Phy -> VIF ; Fun -> VIF ; Beta -> VIF
  VIF -> Brms -> Imp
  Beta -> Maps ; Tax -> Maps ; Phy -> Maps ; Fun -> Maps
  Tax -> Curves ; Phy -> Curves ; Fun -> Curves
  Maps -> Paper ; Curves -> Paper ; Imp -> Paper
}
'

svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
pdf_out <- v3_file("figure", "fig_v3_workflow_pipeline", "pdf")
png_out <- v3_file("figure", "fig_v3_workflow_pipeline", "png")
rsvg::rsvg_pdf(charToRaw(svg), file = pdf_out, width = 1800)
rsvg::rsvg_png(charToRaw(svg), file = png_out, width = 2400)

message("[stage-20] done.")
