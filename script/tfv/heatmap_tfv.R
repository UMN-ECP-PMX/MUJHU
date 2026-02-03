
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(readxl)
library(Tjazi)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/tfv/heatmaps")
if(!dir.exists(figDir)) dir.create(figDir)

theme_set(theme_bw())

# Load and reformat heatmap main panel data -------------------------------

df <- read.csv(here(dataDir, "data.csv"), header = TRUE) %>% 
  mutate(across(contains("phylum"), ~.x/100)) %>% 
  mutate(across(contains("genus"), ~.x/100))

spec <- yspec::load_spec(here(dataDir, "tfv/tfv_expo.yml"))

# Get phylum orders
phylum_order <- df %>% dplyr::select(contains("phylum")) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "phylum", values_to = "mean") %>%
  arrange(desc(mean)) %>% 
  mutate(phylum=gsub("phylum_", "", phylum)) %>% 
  mutate(phylum=str_to_sentence(phylum)) %>% 
  mutate(phylum=fct_inorder(phylum)) %>% 
  pull(phylum) %>% as.factor()

# Select top 20 genus based on abundance
top_genus <- df %>% dplyr::select(contains("genus")) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "genus", values_to = "mean") %>%
  arrange(desc(mean)) %>% 
  slice_head(n=20) %>% pull(genus)

# CLR transformation of relative abundance

# CLR transformation phylum
xx <- df %>% dplyr::select(contains("phylum"))
xx %>% reframe(total=rowSums(across(everything()))) # Make sure sum up to 100 (compositional data)
xx <- clr_c(counts=xx, samples_are = "rows")
xx %>% reframe(total=rowSums(across(everything()))) # Make sure sum up to 0

# CLR transformation genus
yy <- df %>% dplyr::select(contains("genus"))
yy %>% reframe(total=rowSums(across(everything())))
yy <- clr_c(counts=yy, samples_are = "rows")
yy %>% reframe(total=rowSums(across(everything())))

# Do a few plots (CLR transformed phylum/genus vs original)
p1 <- map(names(xx), function(.x){
  p <- ggplot()+geom_point(aes(df[[.x]], xx[[.x]]))+
    xlab(paste0("Original\n", .x))+
    ylab(paste0("CLR transformed\n", .x))
  return(p)
}) %>% pmplots::pm_grid(ncol=6); p1

p2 <- map(names(yy[, top_genus]), function(.x){
  p <- ggplot()+geom_point(aes(df[[.x]], yy[[.x]]))+
    xlab(paste0("Original\n", .x))+
    ylab(paste0("CLR transformed\n", .x))
  return(p)
}) %>% pmplots::pm_grid(ncol=5); p2

# Replace with CLR transformed phylum and genus abundance
df2 <- df %>% dplyr::select(!contains("phylum")) %>%
  dplyr::select(!contains("genus")) %>%
  bind_cols(xx) %>%bind_cols(yy)

# Reformat heatmap data
data <- df2 %>% dplyr::select(contains("tfv"), contains("genus")) %>% 
  dplyr::select(-tfvdp_t, -tfv_p) %>% # Remove raw plasma and cervical concentrations
  dplyr::select(-tfv_ka, -tfv_vc, -tfv_cltt, -tfv_clvvtp, -tfv_clvetp, -tfv_tvkg) %>% # Remove EBEs
  dplyr::select(!(contains("cmax"))) %>% # Remove Cmax (highly correlated with AUC)
  dplyr::select(!(contains("cave"))) # Remove Cmax (highly correlated with AUC)

# Load and reformat top annotation data -----------------------------------

# Extract existing genus in the dataset
exist_genus <- df %>% dplyr::select(contains("genus")) %>% 
  rename_all(.funs = ~gsub("genus_", "", .x)) %>% 
  rename_all(.funs = str_to_sentence) %>% names()

exist_genus_top_20 <- df %>% dplyr::select(contains("genus")) %>% 
  dplyr::select(all_of(top_genus)) %>% # Only filter down to top 20 abundant genus
  rename_all(.funs = ~gsub("genus_", "", .x)) %>% 
  rename_all(.funs = str_to_sentence) %>% names()

# Phylum to Genus alignment for large heatmap annotations
phylum_to_genus_l <- read_xlsx(file.path(here("data/source"), 
                                         "Vaginal Microbiome Data Analysis.xlsx"), 
                      sheet="Taxonomy", col_names=TRUE, range="A1:B494") %>%
  mutate(across(everything(), ~gsub("/", ".", .x))) %>% 
  mutate(across(everything(), str_to_sentence)) %>% 
  rename_all(.funs = str_to_sentence) %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Genus %in% exist_genus)

# Phylum to Genus alignment for small heatmap annotations
phylum_to_genus_s <- phylum_to_genus_l %>% filter(Genus %in% exist_genus_top_20)

# Phylum annotation colors
phylum_colors <- list(Phylum = c(
  "Firmicutes"                  = "#f44336",
  "Actinobacteria"              = "#e81e63",
  "Fusobacteria"                = "#9c27b0",
  "Tenericutes"                 = "#673ab7",
  "Bacteroidetes"               = "#3f51b5",
  "Proteobacteria"              = "#2196f3",
  "Verrucomicrobia"             = "#03a9f4",
  "Candidatus_saccharibacteria" = "#00bcd4",
  # "Sr1"                         = "#009688",
  # "Candidate_division_zb3"      = "#4caf50",
  "Chlamydiae"                  = "#8bc34a",
  # "Cyanobacteria.chloroplast"   = "#cddc39",
  # "Acidobacteria"               = "#ffeb3b",
  # "Candidate_division_wps-1"    = "#ffc107",
  "Planctomycetes"              = "#ff9800",
  "Synergistetes"               = "#ff5722",
  "Elusimicrobia"               = "#4caf50", # "#ffffff",
  "Spirochaetes"                = "#ffeb3b" # "#000000" 
))

# Heatmap function --------------------------------------------------------

hp <- function(dat, top_20_genus=FALSE, phylum_to_genus, d1=TRUE){
  
  dose <- "_d1"
  if (!d1) dose <- "_d2"
  
  # Formatting matrix
  matrix <- cor(dat) %>% as.data.frame() %>% 
    dplyr::select(contains("tfv")) %>% 
    dplyr::select(contains(dose)) %>% 
    rename_all(.fun=~gsub("tfv_", "", .x)) %>% 
    rename_all(.fun=~gsub(dose, "", .x)) %>% 
    rename_all(.fun=toupper) %>% 
    rownames_to_column(var = "names") %>% 
    filter(grepl("genus", names))
  
  if (top_20_genus) matrix <- matrix %>% filter(names %in% top_genus) # Only filter down to top 20 abundant genus
    
  matrix <- matrix %>%  mutate(names = gsub("genus_", "", names)) %>% 
    mutate(names = str_to_sentence(names)) %>% 
    column_to_rownames(var = "names") %>% 
    t() %>% as.matrix()
  
  # Get the desired order of columns (genus) by phylum order
  column_order <- phylum_to_genus %>%
    mutate(Phylum = factor(Phylum, levels = levels(phylum_order))) %>%
    arrange(Phylum) %>% pull(Genus)
  
  # Reorder matrix columns
  matrix <- matrix[, column_order]
  
  # Match the order of phylum (top bar) to the matrix
  phylum_to_genus_out <- phylum_to_genus[match(column_order, phylum_to_genus$Genus), ]
  
  # Output list
  list <- list(matrix, phylum_to_genus_out)
  names(list) <- c("Matrix", "Top-bar")
  
  return(list)
}

# Heat map 2-1-1 dose -----------------------------------------------------

hp1 <- hp(dat=data, top_20_genus=FALSE, 
          phylum_to_genus=phylum_to_genus_l, d1=TRUE)

# Make sure matrix and top-bar genus order aligns
head(colnames(hp1[["Matrix"]]))
head(hp1[["Top-bar"]]$Genus)
tail(colnames(hp1[["Matrix"]]))
tail(hp1[["Top-bar"]]$Genus)

# Heatmap annotations
ha = HeatmapAnnotation(
  Phylum = hp1[["Top-bar"]]$Phylum,
  col = phylum_colors,
  show_annotation_name = c(phylum = FALSE)
)

# Generate heatmap
png(file=here(figDir, "heatmap_d1.png"), width=1800, height=600)
draw(Heatmap(hp1[["Matrix"]], cluster_rows = FALSE, cluster_columns = TRUE, 
             top_annotation = ha, 
             column_split = factor(hp1[["Top-bar"]]$Phylum, levels(phylum_order)), 
             column_order = colnames(hp1[["Matrix"]]), 
             column_title=NULL,
             row_names_gp = grid::gpar(fontsize = 17), 
             column_names_gp = grid::gpar(fontsize = 10), 
             heatmap_legend_param = list(
               legend_direction = "horizontal", 
               title = "")), 
     heatmap_legend_side = "bottom")
dev.off()

# Heat map 2-1-1 dose top20 -----------------------------------------------

hp1s <- hp(dat=data, top_20_genus=TRUE, 
           phylum_to_genus=phylum_to_genus_s, d1=TRUE)

# Make sure matrix and top-bar genus order aligns
head(colnames(hp1s[["Matrix"]]))
head(hp1s[["Top-bar"]]$Genus)
tail(colnames(hp1s[["Matrix"]]))
tail(hp1s[["Top-bar"]]$Genus)

# Heatmap annotations
ha = HeatmapAnnotation(
  Phylum = hp1s[["Top-bar"]]$Phylum,
  col = phylum_colors,
  show_annotation_name = c(phylum = FALSE)
)

# Generate heatmap
png(file=here(figDir, "heatmap_d1s.png"), width=900, height=600)
draw(Heatmap(hp1s[["Matrix"]], cluster_rows = FALSE, cluster_columns = TRUE, 
             top_annotation = ha, 
             column_split = factor(hp1s[["Top-bar"]]$Phylum, levels(phylum_order)), 
             column_order = colnames(hp1s[["Matrix"]]), 
             column_title=NULL,
             row_names_gp = grid::gpar(fontsize = 17), 
             column_names_gp = grid::gpar(fontsize = 17), 
             heatmap_legend_param = list(
               legend_direction = "horizontal", 
               title = "")), 
     heatmap_legend_side = "bottom")
dev.off()

# Heat map 2-1-1-1-1 dose -------------------------------------------------

hp2 <- hp(dat=data, top_20_genus=FALSE, 
          phylum_to_genus=phylum_to_genus_l, d1=FALSE)

# Make sure matrix and top-bar genus order aligns
head(colnames(hp2[["Matrix"]]))
head(hp2[["Top-bar"]]$Genus)
tail(colnames(hp2[["Matrix"]]))
tail(hp2[["Top-bar"]]$Genus)

# Heatmap annotations
ha = HeatmapAnnotation(
  Phylum = hp2[["Top-bar"]]$Phylum,
  col = phylum_colors,
  show_annotation_name = c(phylum = FALSE)
)

# Generate heatmap
png(file=here(figDir, "heatmap_d2.png"), width=1800, height=600)
draw(Heatmap(hp2[["Matrix"]], cluster_rows = FALSE, cluster_columns = TRUE, 
             top_annotation = ha, 
             column_split = factor(hp2[["Top-bar"]]$Phylum, levels(phylum_order)), 
             column_order = colnames(hp2[["Matrix"]]), 
             column_title=NULL,
             row_names_gp = grid::gpar(fontsize = 17), 
             column_names_gp = grid::gpar(fontsize = 10), 
             heatmap_legend_param = list(
               legend_direction = "horizontal", 
               title = "")), 
     heatmap_legend_side = "bottom")
dev.off()

# Heat map 2-1-1-1-1 dose top20 -------------------------------------------

hp2s <- hp(dat=data, top_20_genus=TRUE, 
           phylum_to_genus=phylum_to_genus_s, d1=FALSE)

# Make sure matrix and top-bar genus order aligns
head(colnames(hp2s[["Matrix"]]))
head(hp2s[["Top-bar"]]$Genus)
tail(colnames(hp2s[["Matrix"]]))
tail(hp2s[["Top-bar"]]$Genus)

# Heatmap annotations
ha = HeatmapAnnotation(
  Phylum = hp2s[["Top-bar"]]$Phylum,
  col = phylum_colors,
  show_annotation_name = c(phylum = FALSE)
)

# Generate heatmap
png(file=here(figDir, "heatmap_d2s.png"), width=900, height=600)
draw(Heatmap(hp2s[["Matrix"]], cluster_rows = FALSE, cluster_columns = TRUE, 
             top_annotation = ha, 
             column_split = factor(hp2s[["Top-bar"]]$Phylum, levels(phylum_order)), 
             column_order = colnames(hp2s[["Matrix"]]), 
             column_title=NULL,
             row_names_gp = grid::gpar(fontsize = 17), 
             column_names_gp = grid::gpar(fontsize = 17), 
             heatmap_legend_param = list(
               legend_direction = "horizontal", 
               title = "")), 
     heatmap_legend_side = "bottom")
dev.off()
