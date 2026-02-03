
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

dataDir <- here("data/derived/mn")

figDir <- here("deliv/figure/mn/ftc/pcoa")
if(!dir.exists(figDir)) dir.create(figDir)

theme_set(theme_bw())

options(mrg.script="pcoa_ftc_mn.R", mrggsave.dir=figDir)

# Load absolute count data  -----------------------------------------------

df <- read_xlsx(file.path(here("data/source/mn-biopsy-study"),  
                          "HIV 16S rRNA Data (1).xlsx"), 
                sheet="Genus", col_names=TRUE, range="A1:V265")

# Rename UMN-021 to UMN-022 in the 16s dataset
df <- dplyr::rename(df, MN_022 = MN_021)

# Keep xxxxxaceae_unclassified, remove other _unclassified
df <- df %>%
  filter(!grepl("_unclassified$", Count) | grepl("aceae_unclassified$", Count)) %>%
  column_to_rownames(var = "Count")

# # Initially planned to pre-filter based on percentage of non-zero raw count and abundance 
# # Later found `filter_taxa` function can do it at a later stage
# df <- df %>% rowwise() %>%
#   mutate(
#     fraction_non_zero = sum(c_across(everything()) != 0) / ncol(across(everything()))
#   ) %>%
#   ungroup() %>% 
#   mutate(total = rowSums(.[names(.) != "fraction_non_zero"]))

colnames(df) <- as.numeric(gsub("^MN_", "", colnames(df)))

matrix <- as.matrix(df)

otu <- otu_table(matrix, taxa_are_rows = TRUE)

# Load metadata -----------------------------------------------------------

metadata <- read.csv(here(dataDir, "data.csv"), header = TRUE) %>% 
  dplyr::select(!(contains("genus"))) %>% 
  mutate(ftc_t_exposures = if_else(ftc_auc_t_d1 > median(ftc_auc_t_d1), 
                                   "High", "Low")) %>% 
  mutate(ftc_t_exposures = as.factor(ftc_t_exposures)) %>% 
  mutate(ftc_t_exposures = fct_relevel(ftc_t_exposures, "Low", "High"))

metadata %>% dplyr::count(ftc_t_exposures)

metadata <- metadata %>% column_to_rownames(var = "subject_id")

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
adonis <- adonis2(dist ~ ftc_t_exposures, 
                         data = data.frame(sample_data(ph.obj.rel)))

p_perm <- adonis$`Pr(>F)`[1]
p_label <- paste0("PERMANOVA p = ", signif(p_perm, 3))

# Plot PCoA: Label with ftc exposure
p1 <- plot_ordination(ph.obj.rel, pcoa, color = "ftc_t_exposures") +
  geom_point(size = 3) +
  stat_ellipse() +
  labs(color = "FTC Cervical Exposure",
       subtitle = p_label)
p1
mrggsave(p1, stem="pcoa_plot_genus", width=6, height=5)


# Alpha versus exposure box plots
# Remove id with abnormal chao record  
metadata1 <- metadata %>%
  filter(alpha_chao < 1000)
p1 <- ggboxplot(metadata1, x = "ftc_t_exposures", y = "alpha_shannon", fill = "ftc_t_exposures") +
  stat_compare_means(label.x = 1.5) +
  labs(x = "FTC cervical tissue exposure group", y = "Shannon Diversity")+ 
  theme(legend.position = "none")
p2 <- ggboxplot(metadata1, x = "ftc_t_exposures", y = "alpha_chao", fill = "ftc_t_exposures") +
  stat_compare_means(label.x = 1.5) +
  labs(x = "FTC cervical tissue exposure group", y = "Chao1 Diversity")+ 
  theme(legend.position = "none")
p1 + p2
mrggsave_last(stem="alpha_by_ftc_t_expo", width=8, height=4)
rm(p1,p2)

# Network analysis--------------------------------------------------------

# Extract OTU data by group
#low_otu <- t(as.matrix(otu_table(subset_samples(ph.obj, ftc_t_exposures == "Low"))))
#high_otu <- t(as.matrix(otu_table(subset_samples(ph.obj, ftc_t_exposures == "High"))))

# Rename *_unclassified -> (f)TaxonName
taxa_names(ph.obj) <- taxa_names(ph.obj) %>%
  str_replace("^(.+?)_unclassified$", "(f)\\1")

# Network construction
net_const <- netConstruct(data = ph.obj,
                          #filtTax = "highestVar", filtTaxPar = list(highestVar = 30),
                          measure = "sparcc", #normMethod = "clr", 
                          thresh = 0.3)

# Network analysis
net_props <- netAnalyze(net_const, clustMethod = "cluster_fast_greedy", 
                        hubPar=c("degree","between","closeness"), hubQuant = 0.9)
summary(net_props)

# Plot
pdf(file.path(figDir, "network_plot_genus.pdf"), width = 12, height = 8)
plot(net_props,
     sameLayout = TRUE,
     labelLength = 15,
     labelScale = FALSE,
     nodeFilter="highestEigen",
     nodeFilterPar = 50, 
     rmSingles = FALSE,
     nodeSize="mclr",
     nodeColor = "cluster", 
     hubBorderCol = "black",
     cexNodes = 2,
     cexLabels = 0.8,
     cexHubLabels = 1,
     showTitle = TRUE,
     cexTitle = 1.5)
dev.off()
