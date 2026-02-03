
# Librarys
library(tidyverse)
library(here)
library(mrggsave)
library(readxl)
library(vegan)
library(NetCoMi)
library(ggpubr)
library(ggExtra)
library(patchwork)

# Run these in command line before load `phyloseq`!!
# module load patchelf  # If on a cluster and you need to load it
# patchelf --set-rpath /home/cheng423/shared/software.install/glpk-5.0/lib /panfs/jay/groups/35/cheng423/shared/MUJHU/renv/library/R-4.3/x86_64-pc-linux-gnu/igraph/libs/igraph.so
# After MSI storage migration, the command line should be this:
# patchelf --set-rpath /projects/standard/cheng423/shared/software.install/glpk-5.0/lib /projects/standard/cheng423/shared/MUJHU/renv/library/R-4.3/x86_64-pc-linux-gnu/igraph/libs/igraph.so
library(phyloseq)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/tfv/pcoa")
if(!dir.exists(figDir)) dir.create(figDir)

theme_set(theme_bw())

options(mrg.script="pcoa_tfv.R", mrggsave.dir=figDir)

# Load absolute count data  -----------------------------------------------

df <- read_xlsx(file.path(here("data/source"), 
                          "Vaginal Microbiome Data Analysis.xlsx"), 
                sheet="Genus", col_names=TRUE, range="A1:AY182") %>% 
  column_to_rownames(var = "...1")

matrix <- as.matrix(df)

otu <- otu_table(matrix, taxa_are_rows = TRUE)

# Load metadata -----------------------------------------------------------

metadata <- read.csv(here(dataDir, "data.csv"), header = TRUE) %>% 
  dplyr::select(!(contains("phylum"))) %>% 
  dplyr::select(!(contains("genus"))) %>% 
  mutate(tfv_t_exposures = if_else(tfv_auc_t_d1 > median(tfv_auc_t_d1), 
                                    "High", "Low")) %>% 
  mutate(tfv_t_exposures = as.factor(tfv_t_exposures)) %>% 
  mutate(tfv_t_exposures = fct_relevel(tfv_t_exposures, "Low", "High"))

metadata %>% dplyr::count(tfv_t_exposures)

metadata <- metadata %>% column_to_rownames(var = "msi_id")

samples <- sample_data(metadata)

# Create phyloseq object --------------------------------------------------
ph.obj <- phyloseq(otu, samples)

# PCoA --------------------------------------------------------------------

# Transform to relative abundance
ph.obj.rel <- transform_sample_counts(ph.obj, function(x) x / sum(x))

# Calculate distance matrix using bray curtis
dist <- phyloseq::distance(ph.obj.rel, method = "bray")

# Perform PCoA
pcoa <- ordinate(ph.obj.rel, method = "PCoA", distance = dist)

# Test for significant clustering (PERMANOVA)
adonis <- adonis2(dist ~ tfv_t_exposures, 
                         data = data.frame(sample_data(ph.obj.rel)))

p_perm <- adonis$`Pr(>F)`[1]
p_label <- paste0("PERMANOVA p = ", signif(p_perm, 3))

# Plot PCoA: Label with tfv exposure
p1 <- plot_ordination(ph.obj.rel, pcoa, color = "tfv_t_exposures") +
  geom_point(size = 3) +
  stat_ellipse() +
  labs(color = "TFV Cervical Exposure",
       subtitle = p_label)
p1
mrggsave(p1, stem="pcoa_plot_genus", width=6, height=5)


# Alpha versus exposure box plots
p1 <- ggboxplot(metadata, x = "tfv_t_exposures", y = "alpha_shannon", fill = "tfv_t_exposures") +
  stat_compare_means(label.x = 1.5) +
  labs(x = "TFV cervical tissue exposure group", y = "Shannon Diversity")+ 
  theme(legend.position = "none")
p2 <- ggboxplot(metadata, x = "tfv_t_exposures", y = "alpha_chao", fill = "tfv_t_exposures") +
  stat_compare_means(label.x = 1.5) +
  labs(x = "TFV cervical tissue exposure group", y = "Chao1 Diversity")+ 
  theme(legend.position = "none")
p1 + p2
mrggsave_last(stem="alpha_by_tfv_t_expo", width=8, height=4)
rm(p1,p2)

