library(tidyverse)
library(mrggsave)
library(here)
library(ggvenn)
library(patchwork)

figDir <- here("deliv/figure/mn/ftc/selection")
if(!dir.exists(figDir)) dir.create(figDir)

options(mrggsave.dev="pdf", 
        mrggsave.dir=figDir, 
        mrg.script="selection_ftc_mn.R")

theme_set(theme_bw())

# Load data
df <- tribble(
  ~var,              ~LASSO, ~RF, ~XGB,
  "FTC plasma AUC",       0,     1,    1,
  "Age",                  1,     1,    0,
  "Lachnospiraceae(f)",   0,     1,    1,
  "Gardnerella",          0,     1,    0,
  "Bifidobacteriaceae(f)",1,     1,    1,
  "Prevotellaceae(f)",    0,     1,    1,
  "Mycoplasmopsis",       0,     0,    1,
  "Weight",               0,     0,    1,
  "Prevotella",           1,     0,    0,
  "Megasphaera",          1,     0,    0,
  "Dialister",            1,     0,    0,
  "Mobiluncus",           1,     0,    0
)

lasso_vars <- df$var[df$LASSO == 1]
rf_vars    <- df$var[df$RF    == 1]
xgb_vars   <- df$var[df$XGB   == 1]

set_list <- list(
  LASSO = lasso_vars,
  RF    = rf_vars,
  XGB   = xgb_vars
)

# Heat map

df_long <- df %>%
  mutate(var = factor(var, levels = rev(unique(var)))) %>%
  pivot_longer(
    cols      = c(LASSO, RF, XGB),
    names_to  = "method",
    values_to = "selected"
  )

p <- ggplot(df_long, aes(x = method, y = var, fill = factor(selected))) +
  geom_tile(color = "grey80") +
  scale_fill_manual(
    values = c("0" = "white", "1" = "#2b8cbe"),
    guide  = "none"
  ) +
  scale_x_discrete(position = "top") +
  ggtitle("Selected Variables Across Models") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14))

print(p)

mrggsave_last(
  stem   = "heatmap_selected_genus",
  width  = 5,
  height = 5
)

# Venn plot

p_venn <- ggvenn(
  set_list,
  fill_color      = c("#FFCCBC", "#C5E1A5", "#BBDEFB"),
  stroke_size     = 0.8,
  set_name_size   = 5,
  show_percentage = FALSE
) +
  ggtitle("Overlap of Selected Variables") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  )

print(p_venn)

mrggsave_last(
  stem   = "venn_selected_genus",
  width  = 5,
  height = 5
)

p_fig <- (p + p_venn) +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 14, face = "plain")
  )


print(p_fig)

mrggsave(
  p_fig,
  stem   = "heatmap_venn",
  width  = 10,
  height = 5,
  dpi    = 300,
  compression = "lzw"
)
