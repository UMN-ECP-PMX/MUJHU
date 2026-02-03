
# Librarys
library(tidyverse)
library(here)
library(mrggsave)
library(readxl)
library(ggrepel)
library(DESeq2)

# Run these in command line before load `phyloseq`!!
# module load patchelf  # If on a cluster and you need to load it
# patchelf --set-rpath /home/cheng423/shared/software.install/glpk-5.0/lib /panfs/jay/groups/35/cheng423/shared/MUJHU/renv/library/R-4.3/x86_64-pc-linux-gnu/igraph/libs/igraph.so
library(phyloseq)

dataDir <- here("data/derived/mn")

figDir <- here("deliv/figure/mn/ftc/diffa")
if(!dir.exists(figDir)) dir.create(figDir)

theme_set(theme_bw())

options(mrg.script="diffa_ftc_mn.R", mrggsave.dir=figDir)

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

# Keep taxa present in at least 10% with count of 1 of samples
ph.obj <- filter_taxa(ph.obj, function(x) sum(x > 1) > (0.1 * length(x)), TRUE)

ph.obj
head(sample_data(ph.obj)$ftc_t_exposures, 25)

# Convert a `phyloseq` object to a `deseq2` object ------------------------

deseq.obj <- phyloseq_to_deseq2(ph.obj, ~ ftc_t_exposures)

# Perform DESeq analysis --------------------------------------------------

deseq.obj <- DESeq(deseq.obj, test = "Wald", fitType = "parametric", 
                   sfType = "poscounts") # This handles "0"s in the raw count
                                         # Ignoring "0"s when performing geometric means calculations

# Check results -----------------------------------------------------------

# Gets the results from the object
res <- results(deseq.obj)

# Creates a data frame from results
resdf <- as.data.frame(res) %>% 
  arrange(log2FoldChange, padj) %>%
  # Add labellings for plots
  mutate(sig = ifelse(padj < 0.2 & !is.na(padj), "FDR < 0.2", "FDR >= 0.2")) %>% 
  mutate(sig2 = case_when(
    padj < 0.2 & !is.na(padj) & log2FoldChange < 0 ~ "Negative effect size with FDR < 0.2", 
    padj < 0.2 & !is.na(padj) & log2FoldChange > 0 ~ "Positive effect size with FDR < 0.2", 
    TRUE ~ "FDR >= 0.2"))

# Format genus labels for plotting only ----
fmt <- function(x){
  x <- ifelse(grepl("_unclassified$", x),
              paste0(sub("_unclassified$", "", x), "(f)"),
              x)
  stringr::str_to_sentence(x)
}

rownames(resdf) <- fmt(rownames(resdf))

# Volcano plot ------------------------------------------------------------

# Subset to just the significant ones for labeling
resdf_sig <- resdf %>% filter(sig == "FDR < 0.2")

p1 <- ggplot(resdf, aes(x = log2FoldChange, y = -log10(padj), color = sig2)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Negative effect size with FDR < 0.2" = "blue", 
                                "Positive effect size with FDR < 0.2" = "red", 
                                "FDR >= 0.2" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + # log2FoldChange| > 1 (2-fold change)
  geom_hline(yintercept = -log10(0.2), linetype = "dashed") + # padj < 0.2 (-log10(0.2) ≈ 0.70)
  geom_text_repel(data = resdf_sig, 
                  aes(label = rownames(resdf_sig)), 
                  size = 3, 
                  max.overlaps = 20) +
  labs(x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  #scale_x_continuous(limits=c(-5, 5))+
  theme(legend.position = "none")
p1
mrggsave(p1, stem="volcano_plot_genus", width=6, height=5)

# Scatter plots -----------------------------------------------------------

p2 <- resdf %>% 
  rownames_to_column(var="genus") %>% 
  mutate(genus=fct_inorder(genus)) %>% 
  ggplot(aes(x = genus, y = log2FoldChange, fill = sig2))+
  geom_col(alpha=0.5)+
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Genus", y = "log2 Fold Change")+
  pmplots::rot_x(angle=90)+
  theme(legend.title = element_blank(), legend.position = "bottom")
p2 
mrggsave(p2, stem="scatter_plot_genus", width=15, height=6)
