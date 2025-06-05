#!/usr/bin/env Rscript

# Purpose: What difference would make to use DESeq2 with individual contrast averages

library(tidyverse)
library(ggpubr)

################################################################################
# Load input data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

raw.counts <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

deg30focused <-
  'analysis/K_logFC-vs-30.tsv' |>
  read_tsv()


################################################################################
# Build matrices

raw.counts.mat <-
  raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid)

################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Run DESeq2 with an averaged contrast (version with a model with intercept)

# des <- DESeq2::DESeqDataSetFromMatrix(
#   countData = raw.noribo.mat,
#   colData = meta,
#   design = ~ CO2
# ) |>
#   DESeq2::DESeq()
# 
# 
# deg_check2 <-
#   DESeq2::results(
#     des,
#     contrast = c(
#       # resultsNames(des):
#       # [1] "Intercept"      "CO2_4_vs_0.04"  "CO2_8_vs_0.04"  "CO2_30_vs_0.04"
#       # The contrast vector is derived by:
#       #
#       # 30% condition
#       # Intercept + "CO2_30_vs_0.04"
#       # vs average of 0.04, 4, and 8%
#       #  ((Intercept) + (Intercept + "CO2_4_vs_0.04") + (Intercept + "CO2_8_vs_0.04")) / 3
#       # = ( 3 * Intercept + "CO2_4_vs_0.04" + "CO2_8_vs_0.04") / 3
#       # = Intercept + ("CO2_4_vs_0.04" + "CO2_8_vs_0.04") / 3
#       #
#       # The intercept is in both cases is thus removed, therefore it is
#       "Intercept" = 0,
#       "CO2_4_vs_0.04" = - 1/3,
#       "CO2_8_vs_0.04" = -1/3,
#       "CO2_30_vs_0.04" = 1
#     ),
#     tidy = TRUE
#   ) |>
#   as_tibble() |>
#   rename(Geneid = row)


################################################################################
# Run DESeq2 with an averaged contrast

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  design = ~ CO2 + 0
) |>
  DESeq2::DESeq()


deg_check <-
  DESeq2::results(
    des,
    contrast = c(
      # DESeq2::resultsNames(des)
      # "CO20.04" "CO24"    "CO28"    "CO230" 
      "CO20.04" = - 1/3,
      "CO24" = - 1/3,
      "CO28" = - 1/3,
      "CO230" = + 1
    ),
    tidy = TRUE
  ) |>
  as_tibble() |>
  rename(Geneid = row)

write_tsv(deg_check, 'analysis/N_avg-vs-30.tsv')

################################################################################
# comparison of the base means

cmp_dat <- inner_join(
  deg30focused |>
    select(Geneid, base30 = baseMean),
  deg_check |>
    select(Geneid, baseavg = baseMean),
  'Geneid'
)

cmp_dat |>
  ggscatter('base30', 'baseavg')
# base means are identical

################################################################################
# comparison of the base means

cmp_dat <- inner_join(
  deg30focused |>
    select(Geneid, logFC30 = log2FoldChange, baseMean),
  deg_check |>
    select(Geneid, logFCavg = log2FoldChange),
  'Geneid'
)
  
cmp_dat |>
  mutate_at('baseMean', log10) |>
  ggscatter(
    x = 'logFC30', y = 'logFCavg',
    color = 'baseMean',
    add = 'reg.line',
    add.params = list(color = 'red'),
    cor.coef = TRUE,
    cor.coeff.args = list(color = 'red')
    
  ) +
  scale_color_viridis_c(name = 'log10 average expression (DESeq2 baseMean)') +
  xlab('logFC of 30% CO2 vs all other condition') +
  ylab('logFC of 30% CO2 vs the average of all other conditions') +
  geom_abline(slope = 1, color = 'blue')

ggsave('analysis/N_logFC-scatter-plot.jpeg',
       width = 7, height = 7, dpi = 500)

################################################################################

venn::venn(
  list(
    '30% CO2 vs\nall other condition' =
      deg30focused |>
        filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
        pull(Geneid),
    '30% CO2 vs the\naverage of all other conditions' =
      deg_check |>
        filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
        pull(Geneid)
  ),
  ilabels = 'counts',
  zcolor = 'style',
  ilcs = 1,
  sncs = 1,
  ggplot = TRUE,
  box = FALSE
) +
  annotate(
   'text', 500, 900,
   label = 'DEG for abs logFC ≥ 1 and FDR 5%'
  )

ggsave('analysis/N_venn.jpeg',
       width = 7, height = 7, dpi = 500)

################################################################################

xs <- symdiff(
  deg30focused |>
    filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
    pull(Geneid),
  deg_check |>
    filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
    pull(Geneid)
)
x_focused <-
  deg30focused |>
    filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
    pull(Geneid)

annot |>
  filter(Geneid %in% xs) |>
  mutate(
    which = ifelse(
      Geneid %in% x_focused,
      '30% vs all',
      '30$ vs avg'
  )) |>
  View()
