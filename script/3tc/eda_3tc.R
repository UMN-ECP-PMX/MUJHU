
library(tidyverse)
library(here)
library(patchwork)
library(mrggsave)
library(ggpubr)
library(ggExtra)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/3tc/eda")
if(!dir.exists(figDir)) dir.create(figDir)

options(mrggsave.dev="pdf", 
        mrggsave.dir=figDir, 
        mrg.script="eda_3tc.R")

theme_set(theme_bw())

# Load data ---------------------------------------------------------------

df <- read.csv(here(dataDir, "data.csv"), header = TRUE)

# Tissue exposure versus contraceptives box plots -------------------------

box_by <- function(var1, var2){
  p <- ggboxplot(df, x=var1, y=var2)+
    stat_compare_means(label.x = 1.5)
  return(p)}

plot_box <- function(dd,xx){
  p      <- map(names(dplyr::select(df, contains("3tc")&contains(paste0("p_", dd)))), 
                ~box_by(var1=xx, var2=.x))
  t      <- map(names(dplyr::select(df, contains("3tc")&contains(paste0("t_", dd)))), 
                ~box_by(var1=xx, var2=.x))

  p      <- pmplots::pm_grid(p, ncol=3)
  t      <- pmplots::pm_grid(t, ncol=3)

  p <- p / t 
  return(p)
}

plot_box(dd="d1",xx="hormone")
mrggsave_last(stem="expo_by_hormone_d1", width=12, height=8)
plot_box(dd="d2",xx="hormone")
mrggsave_last(stem="expo_by_hormone_d2", width=12, height=8)
plot_box(dd="d1",xx="contra")
mrggsave_last(stem="expo_by_contra_d1", width=12, height=8)
plot_box(dd="d2",xx="contra")
mrggsave_last(stem="expo_by_contra_d2", width=12, height=8)

# Tissue exposure versus other categorical box plots ----------------------

df %>% count(sexually_active)
p1 <- plot_box(dd="d1",xx="sexually_active")
p2 <- plot_box(dd="d2",xx="sexually_active")
df %>% count(ammenorheic)
p3 <- plot_box(dd="d1",xx="ammenorheic")
p4 <- plot_box(dd="d2",xx="ammenorheic")
df %>% count(urinary_gonorrhea)
p5 <- plot_box(dd="d1",xx="urinary_gonorrhea")
p6 <- plot_box(dd="d2",xx="urinary_gonorrhea")
df %>% count(rpr_results)
p7 <- plot_box(dd="d1",xx="rpr_results")
p8 <- plot_box(dd="d2",xx="rpr_results")

mrggsave(list(p1,p3,p5,p7), stem="expo_by_cat_d1", width=12, height=8)
mrggsave(list(p2,p4,p6,p8), stem="expo_by_cat_d2", width=12, height=8)

rm(p1,p2,p3,p4,p5,p6,p7,p8)

# Exposure versus alpha box plots -----------------------------------------

df2 <- df %>% mutate(alpha_shannon2 = ifelse(alpha_shannon > median(df$alpha_shannon), 
                                             "High", "Low"), 
                     alpha_chao2 = ifelse(alpha_chao > median(df$alpha_chao), 
                                          "High", "Low"))

ggplot(df2, aes(x=alpha_shannon2, y=X3tc_auc_t_d1))+
  geom_boxplot()+
  facet_wrap(~hormone)

ggplot(df2, aes(x=alpha_chao2, y=X3tc_auc_t_d1))+
  geom_boxplot()+
  facet_wrap(~hormone)

rm(df2)

# Tissue exposure versus alpha scatter plot -------------------------------

plot_scatter <- function(.x, .y, .xlab, .ylab, .color){
  p <- ggplot(df, aes(x={{.x}}, y={{.y}}, 
                      group={{.color}}, color={{.color}}, fill={{.color}}))+
    geom_point()+geom_smooth(method="lm")+
    xlab(.xlab)+ylab(.ylab)+
    theme(legend.position = "bottom")
  p <- ggMarginal(p, type = "boxplot", groupFill = TRUE); 
  return(p)
}

# Stratified by `hormone`
p1 <- plot_scatter(.x=alpha_shannon, .y=X3tc_auc_t_d1, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p1
p2 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cmax_t_d1, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p2
p3 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cave_t_d1, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p3
p4 <- plot_scatter(.x=alpha_chao, .y=X3tc_auc_t_d1, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p4
p5 <- plot_scatter(.x=alpha_chao, .y=X3tc_cmax_t_d1, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p5
p6 <- plot_scatter(.x=alpha_chao, .y=X3tc_cave_t_d1, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p6

mrggsave(list(p1,p2,p3,p4,p5,p6), 
         stem="expo_vs_alpha_by_hormone_d1", 
         width=6, height=6)

p7 <- plot_scatter(.x=alpha_shannon, .y=X3tc_auc_t_d2, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p7
p8 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cmax_t_d2, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p8
p9 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cave_t_d2, .color=hormone,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p9
p10<- plot_scatter(.x=alpha_chao, .y=X3tc_auc_t_d2, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p10
p11<- plot_scatter(.x=alpha_chao, .y=X3tc_cmax_t_d2, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p11
p12<- plot_scatter(.x=alpha_chao, .y=X3tc_cave_t_d2, .color=hormone,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p12

mrggsave(list(p7,p8,p9,p10,p11,p12), 
         stem="expo_vs_alpha_by_hormone_d2", 
         width=6, height=6)

rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)

# Stratified by `contra`
p1 <- plot_scatter(.x=alpha_shannon, .y=X3tc_auc_t_d1, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p1
p2 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cmax_t_d1, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p2
p3 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cave_t_d1, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p3
p4 <- plot_scatter(.x=alpha_chao, .y=X3tc_auc_t_d1, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p4
p5 <- plot_scatter(.x=alpha_chao, .y=X3tc_cmax_t_d1, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p5
p6 <- plot_scatter(.x=alpha_chao, .y=X3tc_cave_t_d1, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p6

mrggsave(list(p1,p2,p3,p4,p5,p6), 
         stem="expo_vs_alpha_by_contra_d1", 
         width=6, height=6)

p7 <- plot_scatter(.x=alpha_shannon, .y=X3tc_auc_t_d2, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p7
p8 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cmax_t_d2, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p8
p9 <- plot_scatter(.x=alpha_shannon, .y=X3tc_cave_t_d2, .color=contra,
                   .xlab="Shannon Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p9
p10<- plot_scatter(.x=alpha_chao, .y=X3tc_auc_t_d2, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue AUC (fmol*hour/mL)"); p10
p11<- plot_scatter(.x=alpha_chao, .y=X3tc_cmax_t_d2, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cmax (fmol/mL)"); p11
p12<- plot_scatter(.x=alpha_chao, .y=X3tc_cave_t_d2, .color=contra,
                   .xlab="Chao Index", 
                   .ylab="3TC viginal tissue Cave (fmol/mL)"); p12

mrggsave(list(p7,p8,p9,p10,p11,p12), 
         stem="expo_vs_alpha_by_contra_d2", 
         width=6, height=6)

rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)


# Microbiome composition plots (phylum level) -----------------------------

# Data processing
data <- df %>% 
  dplyr::select(subject_id, age, weight, hormone, contra, contains("alpha"), 
                contains("X3tc_auc_t"), contains("X3tc_cmax_t"), contains("X3tc_cave_t"), 
                contains("phylum")) %>% 
  pivot_longer(cols=contains("phylum"), names_to="Phylum", values_to="Abundance") %>% 
  mutate(`Phylum`=gsub("phylum_", "", `Phylum`)) %>% 
  mutate(`Phylum`=str_to_sentence(`Phylum`))

p1 <- data %>% 
  ggplot(aes(x=as.factor(X3tc_auc_t_d1), y=Abundance, group=Phylum, fill=Phylum))+
  geom_col(position="stack")+pmplots::rot_x(angle = 60)+
  xlab("3TC cervical tissue exposure low to high")+
  theme(legend.position = "bottom", axis.text.x=element_blank())

p2 <- data %>% 
  ggplot(aes(x=as.factor(alpha_shannon), y=Abundance, group=Phylum, fill=Phylum))+
  geom_col(position="stack")+pmplots::rot_x(angle = 60)+
  xlab("Shannon index low to high")+
  theme(legend.position = "bottom", axis.text.x=element_blank())

mrggsave(list(p1, p2), 
         stem="microbiome_composition_phylum", 
         width=15, height=8)

rm(data, p1, p2)

# Microbiome composition plots (genus level) ------------------------------

# Select top 20 genus based on abundance
top_genus <- df %>% dplyr::select(contains("genus")) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "genus", values_to = "mean") %>%
  arrange(desc(mean)) %>% 
  slice_head(n=20) %>% 
  pull(genus)

# Data processing
data <- df %>% 
  dplyr::select(subject_id, age, weight, hormone, contra, contains("alpha"), 
                contains("X3tc_auc_t"), contains("X3tc_cmax_t"), contains("X3tc_cave_t"), 
                contains("genus")) %>% 
  pivot_longer(cols=contains("genus"), names_to="Genus", values_to="Abundance") %>% 
  filter(Genus %in% top_genus) %>% # Only keep top 20 genus in terms of abundance
  mutate(`Genus`=gsub("genus_", "", `Genus`)) %>% 
  mutate(`Genus`=str_to_sentence(`Genus`)) 

other <- data %>% # Abundance ranks >20 aggreagtes in a single "Other" genera
  group_by(subject_id) %>% 
  summarise(Genus="Other", Abundance=100-sum(Abundance)) %>% 
  ungroup()

data <- full_join(data, other) %>% # Combine two parts together
  group_by(subject_id) %>% 
  fill(age,weight,hormone,contra,alpha_shannon,alpha_chao,
       X3tc_auc_t_d1,X3tc_auc_t_d2,X3tc_cmax_t_d1,X3tc_cmax_t_d2,X3tc_cave_t_d1,X3tc_cave_t_d2,
       .direction = "downup") %>% 
  ungroup()
  
p1 <- data %>% 
  ggplot(aes(x=as.factor(X3tc_auc_t_d1), y=Abundance, group=Genus, fill=Genus))+
  geom_col(position="stack")+pmplots::rot_x(angle = 60)+
  xlab("3TC cervical tissue exposure low to high")+
  theme(legend.position = "bottom", axis.text.x=element_blank())

p2 <- data %>% 
  ggplot(aes(x=as.factor(alpha_shannon), y=Abundance, group=Genus, fill=Genus))+
  geom_col(position="stack")+pmplots::rot_x(angle = 60)+
  xlab("Shannon index low to high")+
  theme(legend.position = "bottom", axis.text.x=element_blank())

mrggsave(list(p1, p2), 
         stem="microbiome_composition_genus", 
         width=15, height=8)


col_every_genus <- function(genus, by){
  
  data2 <- data %>% mutate(Genus=ifelse(Genus==genus, genus, "Everything else")) %>% 
    mutate(Genus = fct_relevel(Genus, c("Everything else", genus)))
  
  if (by == "exposures"){
    p <- data2 %>% ggplot(aes(x=as.factor(X3tc_auc_t_d1)))+
      xlab("3TC cervical tissue exposure low to high")
  } else if (by == "alpha"){
    p <- data2 %>% ggplot(aes(x=as.factor(alpha_shannon)))+
      xlab("Shannon index low to high")
  } else {
    print("`by` can only be `exposures` or `alpha`")
  }
  
  p <- p + geom_col(aes(y=Abundance, group=Genus, fill=Genus), position="stack")+
    pmplots::rot_x(angle = 60)+
    theme(legend.position = "bottom", axis.text.x=element_blank())+
    scale_fill_manual(values = c("red", "gray"), 
                      breaks = c(genus, "Everything else"))  
  
  return(p)
}

p3list <- map(unique(data$Genus), ~col_every_genus(.x, by = "exposures"))
mrggsave(p3list, stem="microbiome_composition_genus_expo", width=15, height=8)

p4list <- map(unique(data$Genus), ~col_every_genus(.x, by = "alpha"))
mrggsave(p4list, stem="microbiome_composition_genus_alpha", width=15, height=8)

