
# Library and settings ----------------------------------------------------

# Load library
library(here)
library(tidyverse)
library(mrgsolve)
library(mapbayr)
library(mrggsave)
library(forcats)

outDatDir <- here("data/derived/mn/ftc")
if (!dir.exists(outDatDir)) dir.create(outDatDir)
outFigDir <- here("deliv/figure/mn/ftc")
if (!dir.exists(outFigDir)) dir.create(outFigDir)

options(mrgsolve.project=here("model/mrgsolve"), 
        mrggsave.dir = outFigDir, 
        mrggsave.dev = "pdf", 
        mrg.script = "mapbayes_ftc_mn.R")

theme_set(theme_bw())

# Data assembly -----------------------------------------------------------

# Load observations
obs <- read.csv(here("data/derived/mn/pkdat.csv")) %>% 
  mutate(blq = ifelse(blq==1, "BLQ", "Non-BLQ"))

# Check observations
obs %>% count(blq)
obs %>% count(flag)
obs %>% count(flag, blq)

obs %>% filter(blq == "BLQ") %>% 
  filter(flag %in% c("ftc_p", "ftctp_t")) %>% 
  count(subject_id, flag)

# Data assembly for observations
obs1 <- obs %>% 
  filter(flag %in% c("ftc_p", "ftctp_t")) %>% 
  # filter(blq != "BLQ") %>% # decided to include BLQ Concs
  filter(!is.na(dv)) %>% 
  mutate(evid=0, amt=0, ss=0, ii=0) %>% 
  mutate(cmt=case_when(
    flag == "ftc_p"~2, 
    flag == "ftctp_t"~8, 
    TRUE ~ NA_real_)) %>%
  mutate(mdv = 0) %>% 
  mutate(ID = as.numeric(as.factor(subject_id))) %>% 
  dplyr::select(ID, time, evid, cmt, amt, DV=dv, mdv, ss, ii, LLOQ=lloq, 
                subject_id, flag)

# Align observation units
obs_final <- obs1 %>% 
  mutate(DV = case_when(
    flag == "ftc_p"~((DV/1e6)/247.25)*1e15, 
    #' ng/ml to fmol/L (ng/mL to ug/L,    # /1
    #'                  ug/L to g/L,      # /1e6
    #'                  g/L to mol/L,     # /247.25 g/mol (ref:https://www.medchemexpress.com/emtricitabine-standard.html)
    #'                  mol/L to fmol/L)  # *1e15
    flag == "ftctp_t"~DV*1e3, # fmol/g=fmol/mL to fmol/L
    TRUE ~ NA_real_)) %>% 
  mutate(LLOQ = case_when(
    flag == "ftc_p"~((LLOQ/1e6)/247.25)*1e15, 
    flag == "ftctp_t"~LLOQ*1e3, # fmol/mL to fmol/L
    TRUE ~ NA_real_)) %>% 
  dplyr::select(-subject_id, -flag) # Remove character cols

# Create doses
dose <- data.frame(
  ID = unique(obs_final$ID), 
  time=0, evid=1, cmt=1, 
  amt = ((200/1000)/247.25)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
  DV=0, mdv=1, ss=1, ii=24, LLOQ=0
)

dose_final <- dose

# Combine dose and observations
data <- rbind(obs_final, dose_final) %>% 
  arrange(ID, time, desc(evid), cmt) 

# Load mrgsolve model -----------------------------------------------------

mod <- mread("3tc_mapbayr.cpp")
mod <- update(mod, atol=1e-60, delta=0.1)

mod
param(mod)
revar(mod)

# Map bayes VPC -----------------------------------------------------------

# Create a labeller for facet renaming
rename_facet <- c(
  PAR = "FTC plasma",
  MET = "FTCtp cervical tissue"
)

# Extract BLQs
obs_all <- data %>% filter(evid == 0) %>%
  mutate(name = case_when(cmt == 2 ~ "PAR", cmt == 8 ~ "MET")) %>% 
           mutate(name = as.factor(name))

blq <- data %>% filter(DV < LLOQ) %>% 
  mutate(name = case_when(cmt==2~"PAR", cmt==8~"MET")) %>% 
  mutate(name = as.factor(name))

# Create VPCs
withr::with_seed(seed = 123, 
                 vpc <- mapbayr::mapbayr_vpc(
                   mod, data = data, nrep = 200, 
                   pcvpc = FALSE, end = 48, delta = 0.1)
)

# Delete obs points created by VPCs
point_idx <- vapply(
  vpc$layers,
  function(ly) inherits(ly$geom, "GeomPoint"),
  logical(1)
)
vpc$layers <- vpc$layers[!point_idx]

# Add obs and blq points on VPC
vpc <- vpc +
  geom_point(data = obs_all, aes(x = time, y = DV)) +
  geom_point(data = blq, aes(x = time, y = DV), color = "red") +
  facet_wrap(~ name, scales  = "free",
    labeller = labeller(name = rename_facet)) +
  xlab("Time after last dose (hour)") +
  ylab("Concentration (fmol/L)")

mrggsave(vpc,stem = "vpc", width=6, height=5)

# Create VPCs in log scale
vpc_log <- vpc + scale_y_log10()
mrggsave(vpc_log,stem = "vpc_log", width=6, height=5)


# Map bayes estimates for all ---------------------------------------------

my_est <- mapbayest(mod, data = data)

# Check estimates
print(my_est)

# Individual plots

mapbayr_plot2 <- function(i, logy = FALSE) {
  
  my_est_tab <- augment(my_est)
  
  aug_tab1 <- my_est_tab$aug_tab %>% 
    filter(ID == i)
  
  obs_tab1 <- my_est_tab$mapbay_tab %>% 
    filter(ID == i, cmt %in% c(2, 8)) %>% 
    mutate(name = case_when(cmt == 2 ~ "PAR",
                            cmt == 8 ~ "MET"),
           name = factor(name, levels = c("PAR", "MET"))
    )
  
  rename_facet <- c(
    PAR = "FTC plasma",
    MET = "FTCtp cervical tissue"
  )
  
  subj_id <- filter(obs1, ID == i) %>% pull(subject_id) %>% unique()
  
  # Only plot PRED and IPRED, without obs
  empty_obs <- obs_tab1[0, ]
  
  p <- mapbayr::mapbayr_plot(aug_tab = aug_tab1, 
                             obs_tab = empty_obs) +
    facet_wrap(~name, scales = "free", 
               labeller = labeller(name=rename_facet))+
    ggtitle(paste0("Subject ID: ", subj_id))+
    xlab("Time after last dose (hour)")+
    ylab("Concentration (fmol/L)")
  
  # Add obs points
  p <- p +
    ggplot2::geom_point(
      data = obs_tab1,
      aes(x = time, y = DV))
  
  ## Add blq points
  blq_i <- obs_tab1 %>% dplyr::filter(DV < LLOQ)
  
  p <- p +
    ggplot2::geom_point(
      data = blq_i,
      aes(x = time, y = DV),
      colour = "red"
    )
  
  if (logy) p <- p + scale_y_log10() 
  
  return(p)
  
}

ind_plots <- map(unique(data$ID), mapbayr_plot2)
mrggsave(ind_plots,stem = "individual_fits", width=6, height=5)

ind_plots_log <- map(unique(data$ID), ~mapbayr_plot2(.x, logy=TRUE))
mrggsave(ind_plots_log,stem = "individual_fits_log", width=6, height=5)

# Histogram
hist <- hist(my_est)
mrggsave(hist,stem = "histogram", width=6, height=5)

# Output Map Bayes estimates ----------------------------------------------

# Extract subject ID
subj_id <- obs1 %>% distinct(ID, subject_id) %>% rename_with(toupper)

# Extract subject ID with BLQs
ids <- obs %>% filter(blq == "BLQ") %>% 
  filter(flag %in% c("ftc_p", "ftctp_t")) %>% 
  count(subject_id, flag) %>% 
  distinct(subject_id) %>% pull(subject_id)

# Extract ETAs
etas <- augment(my_est)$mapbay_tab %>% 
  dplyr::select(ID, contains("ETA")) %>% 
  distinct() %>% 
  # Append `subject_id` in the original data
  left_join(subj_id) %>% 
  # Give a flag for potential exclusion
  mutate(EX_FLAG = if_else(SUBJECT_ID %in% ids, 1, 0)) %>% 
  dplyr::select(SUBJECT_ID, ID, EX_FLAG, everything())

write.csv(etas, file.path(outDatDir, "ftc_etas.csv"), 
          row.names = FALSE, quote = FALSE)
