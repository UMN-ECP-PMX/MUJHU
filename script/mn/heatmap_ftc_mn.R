
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(readxl)
library(Tjazi)
library(viridisLite)

dataDir <- here("data/derived/mn")

figDir <- here("deliv/figure/mn/ftc/heatmaps")
if(!dir.exists(figDir)) dir.create(figDir)

theme_set(theme_bw())

# Load and reformat heatmap main panel data -------------------------------

df <- read.csv(here(dataDir, "data.csv"), header = TRUE) %>% 
  mutate(across(contains("genus"), ~.x/100))

spec <- yspec::load_spec(here(dataDir, "ftc/ftc_expo.yml"))

# Select top 20 genus based on abundance
top_genus <- df %>% dplyr::select(contains("genus")) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "genus", values_to = "mean") %>%
  arrange(desc(mean)) %>% 
  slice_head(n=20) %>% pull(genus)

# CLR transformation of relative abundance

# CLR transformation genus
yy <- df %>% dplyr::select(contains("genus"))
yy %>% reframe(total=rowSums(across(everything())))
yy <- clr_c(counts=yy, samples_are = "rows")
yy %>% reframe(total=rowSums(across(everything())))

# Do a few plots (CLR transformed genus vs original)
p2 <- map(names(yy[, top_genus]), function(.x){
  p <- ggplot()+geom_point(aes(df[[.x]], yy[[.x]]))+
    xlab(paste0("Original\n", .x))+
    ylab(paste0("CLR transformed\n", .x))
  return(p)
}) %>% pmplots::pm_grid(ncol=5); p2

# Replace with CLR transformed genus abundance
df2 <- df %>% 
  dplyr::select(!contains("genus")) %>%
  bind_cols(yy)

# Reformat heatmap data
data <- df2 %>% dplyr::select(contains("ftc"), contains("genus")) %>% 
  dplyr::select(-ftctp.fmol.g : -ftc.plasma) %>% # Remove raw plasma and cervical concentrations
  dplyr::select(-ftc_ka : -ftc_kg) %>% # Remove EBEs
  dplyr::select(!(contains("cmax"))) %>% # Remove Cmax (highly correlated with AUC)
  dplyr::select(!(contains("cave"))) %>% # Remove Cmax (highly correlated with AUC)
  dplyr::select(!(contains("frac24"))) # Remove fraciton above thereshold (all are 1)

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
phylum_to_genus_l <- read_xlsx(file.path(here("data/source/mn-biopsy-study"), 
                                         "HIV 16S rRNA Data (1).xlsx"), 
                      sheet="Taxonomy (norm=30000)", col_names=TRUE, range="A1:AB3800") %>%
  select(Phylum, Genus) %>%
  mutate(across(everything(), ~gsub("[/-]", ".", .x))) %>% 
  mutate(across(everything(), str_to_sentence)) %>% 
  rename_all(.funs = str_to_sentence) %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Genus %in% exist_genus)

# Phylum to Genus alignment for small heatmap annotations
phylum_to_genus_s <- phylum_to_genus_l %>% filter(Genus %in% exist_genus_top_20)
unique(phylum_to_genus_l$Phylum)

# Phylum annotation colors
phylum_list <- c(
  "Firmicutes","Actinobacteria","Fusobacteria","Bacteroidetes","Tenericutes",
  "Campilobacterota","Proteobacteria","Bacteria_unclassified","Archaea_unclassified",
  "Cyanobacteria.chloroplast","Euryarchaeota","Crenarchaeota","Verrucomicrobia",
  "Nitrospirae","Aminicenantes","Unknown_unclassified","Chloroflexi",
  "Woesearchaeota","Planctomycetes","Acidobacteria","Pacearchaeota",
  "Aenigmarchaeota","Candidate_division_wps.1","Candidatus_saccharibacteria",
  "Thaumarchaeota","Armatimonadetes","Latescibacteria","Spirochaetes",
  "Fibrobacteres","Brc1","Synergistetes","Deferribacteres","Hydrogenedentes",
  "Deinococcus.thermus"
)

colors <- turbo(length(phylum_list))

phylum_colors <- list(
  Phylum = setNames(colors, phylum_list)
)
phylum_colors

phylum_order <- factor(phylum_list, levels = phylum_list)

# Heatmap function --------------------------------------------------------

hp <- function(dat, top_20_genus=FALSE, phylum_to_genus, d1=TRUE){
  
  dose <- "_d1"
  if (!d1) dose <- "_d2"
  
  # unclassified to (f) 
  fmt <- function(x){
    x <- gsub("^genus_", "", x)
    x <- ifelse(grepl("_unclassified$", x),
                paste0(sub("_unclassified$", "", x), "(f)"),
                x)
    stringr::str_to_sentence(x)
  }
  
  # Formatting matrix
  matrix <- cor(dat) %>% as.data.frame() %>% 
    dplyr::select(contains("ftc")) %>% 
    dplyr::select(contains(dose)) %>% 
    rename_all(.fun=~gsub("ftc_", "", .x)) %>% 
    rename_all(.fun=~gsub(dose, "", .x)) %>% 
    rename_all(.fun=toupper) %>% 
    rownames_to_column(var = "names") %>% 
    filter(grepl("genus", names))
  
  if (top_20_genus) matrix <- matrix %>% filter(names %in% top_genus) # Only filter down to top 20 abundant genus
    
  matrix <- matrix %>%  
    mutate(names = fmt(names)) %>% 
    column_to_rownames(var = "names") %>% 
    t() %>% as.matrix()
  
  # Get the desired order of columns (genus) by phylum order
  column_order <- phylum_to_genus %>%
    mutate(Phylum = factor(Phylum, levels = levels(phylum_order))) %>%
    arrange(Phylum) %>% 
    mutate(Genus = fmt(Genus)) %>% 
    pull(Genus)
  
  # Reorder matrix columns
  matrix <- matrix[, column_order]
  
  # Match the order of phylum (top bar) to the matrix
  phylum_to_genus_out <- phylum_to_genus %>%
    mutate(Genus = fmt(Genus)) %>%
    .[match(column_order, .$Genus), ]
  
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
png(file=here(figDir, "heatmap_d1.png"), width=2400, height=600)
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
             column_names_gp = grid::gpar(fontsize = 14), 
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
png(file=here(figDir, "heatmap_d2.png"), width=2400, height=600)
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
             column_names_gp = grid::gpar(fontsize = 14), 
             heatmap_legend_param = list(
               legend_direction = "horizontal", 
               title = "")), 
     heatmap_legend_side = "bottom")
dev.off()
