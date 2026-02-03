
# Library and settings ----------------------------------------------------

# Load library
library(here)
library(tidyverse)
library(mrgsolve)
library(mapbayr)
library(mrggsave)

outDatDir <- here("data/derived/3tc")
if (!dir.exists(outDatDir)) dir.create(outDatDir)
outFigDir <- here("deliv/figure/3tc")
if (!dir.exists(outFigDir)) dir.create(outFigDir)

options(mrgsolve.project=here("model/mrgsolve"), 
        mrggsave.dir = outFigDir, 
        mrggsave.dev = "pdf", 
        mrg.script = "mapbayes_3tc.R")

theme_set(theme_bw())

# Data assembly -----------------------------------------------------------

# Load observations
obs <- read.csv(here("data/derived/pkdat.csv")) %>% 
  mutate(blq = ifelse(blq==1, "BLQ", "Non-BLQ"))

# Check observations
obs %>% count(blq)
obs %>% count(flag)
obs %>% count(flag, blq)

obs %>% filter(blq == "BLQ") %>% 
  filter(flag %in% c("_3tc_p", "_3tctp_t")) %>% 
  count(subject_id, flag)

# Data assembly for observations
obs1 <- obs %>% 
  filter(flag %in% c("_3tc_p", "_3tctp_t")) %>% 
  # filter(blq != "BLQ") %>% # decided to include BLQ Concs
  filter(!is.na(dv)) %>% 
  mutate(evid=0, amt=0, ss=0, ii=0) %>% 
  mutate(cmt=case_when(
    flag == "_3tc_p"~2, 
    flag == "_3tctp_t"~8, 
    TRUE ~ NA_real_)) %>%
  mutate(mdv = 0) %>% 
  mutate(ID = as.numeric(as.factor(subject_id))) %>% 
  dplyr::select(ID, time, evid, cmt, amt, DV=dv, mdv, ss, ii, LLOQ=lloq, 
                subject_id, flag)

# Align observation units
obs_final <- obs1 %>% 
  mutate(DV = case_when(
    flag == "_3tc_p"~((DV/1e6)/229.26)*1e15, 
    #' ng/ml to fmol/L (ng/mL to ug/L,    # /1
    #'                  ug/L to g/L,      # /1e6
    #'                  g/L to mol/L,     # /229.26 g/mol (ref:https://www.medchemexpress.com/Lamivudine.html)
    #'                  mol/L to fmol/L)  # *1e15
    flag == "_3tctp_t"~DV*1e3, # fmol/g=fmol/mL to fmol/L
    TRUE ~ NA_real_)) %>% 
  mutate(LLOQ = case_when(
    flag == "_3tc_p"~((LLOQ/1e6)/229.26)*1e15, 
    flag == "_3tctp_t"~LLOQ*1e3, # fmol/mL to fmol/L
    TRUE ~ NA_real_)) %>% 
  dplyr::select(-subject_id, -flag) # Remove character cols

# Create doses
dose <- data.frame(
  ID = unique(obs_final$ID), 
  time=0, evid=1, cmt=1, 
  amt = ((300/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
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
  PAR = "3TC plasma",
  MET = "3TCtp cervical tissue"
)

# Extract BLQs
blq <- data %>% filter(DV < LLOQ) %>% 
  mutate(name = case_when(cmt==2~"PAR", cmt==8~"MET")) %>% 
  mutate(name = as.factor(name))

# Create VPCs
withr::with_seed(seed = 123, 
                 vpc <- mapbayr::mapbayr_vpc(
                   mod, data = data, nrep = 200, 
                   pcvpc = FALSE, end = 48, delta = 0.1)
                 )

vpc <- vpc + geom_point(data=blq, aes(x=time, y=DV), color="red")+
  facet_wrap(~name, scales = "free", 
             labeller = labeller(name=rename_facet))+
  xlab("Time after last dose (hour)")+
  ylab("Concentration (fmol/L)")

mrggsave(vpc,stem = "vpc", width=6, height=5)

# Create VPCs in log scale
vpc_log <- vpc + scale_y_log10()
mrggsave(vpc_log,stem = "vpc_log", width=6, height=5)

# Map bayes individual test -----------------------------------------------

# Test one individual
test_est <- mod %>% 
  adm_rows(ID=1, time=0, cmt=1, amt=472054380665, ss=1, ii=24) %>% 
  obs_rows(ID=1, time=19.666667, cmt=2, DV=208553841) %>% # QC: Looks like TFV values
  obs_rows(ID=1, time=19.933333, cmt=8, DV=9866000) %>% 
  mapbayest()

# Check fittings
print(test_est)
plot(test_est)
hist(test_est)

# Map bayes estimates for all ---------------------------------------------

my_est <- mapbayest(mod, data = data)

# QC: A few warning messages? 
# Reset with new bounds (lower displayed): -1.81006 -0.754566 -1.37593 -1.42943 -1.24374 -4.26231 -2.79187 -2.35128
# Reset with new bounds (lower displayed): -2.07574 -0.86532 -1.57789 -1.63924 -1.4263 -4.88792 -3.20165 -2.6964
# Reset with new bounds (lower displayed): -2.31351 -0.964441 -1.75863 -1.82702 -1.58968 -5.44782 -3.5684 -3.00526
# Reset with new bounds (lower displayed): -2.53053 -1.05491 -1.9236 -1.99841 -1.7388 -5.95888 -3.90314 -3.28718
# Reset with new bounds (lower displayed): -2.73138 -1.13864 -2.07628 -2.15702 -1.87681 -6.43183 -4.21293 -3.54808
# Reset with new bounds (lower displayed): -2.91915 -1.21692 -2.21902 -2.3053 -2.00583 -6.87399 -4.50256 -3.792

# Check estimates
print(my_est)

# Individual plots
mapbayr_plot2 <- function(i, logy=FALSE){
  
  my_est_tab <- augment(my_est)
  
  aug_tab1 <- my_est_tab$aug_tab %>% 
    filter(ID == i)
  
  obs_tab1 <- my_est_tab$mapbay_tab %>% 
    filter(ID == i) 
  
  rename_facet <- c(
    PAR = "3TC plasma",
    MET = "3TCtp cervical tissue"
  )
  
  subj_id <- filter(obs1, ID == i) %>% pull(subject_id) %>% unique()
  
  p <- mapbayr::mapbayr_plot(aug_tab = aug_tab1, 
                        obs_tab = obs_tab1)+
    facet_wrap(~name, scales = "free", 
               labeller = labeller(name=rename_facet))+
    ggtitle(paste0("Subject ID: ", subj_id))+
    xlab("Time after last dose (hour)")+
    ylab("Concentration (fmol/L)")
  
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
  filter(flag %in% c("_3tc_p", "_3tctp_t")) %>% 
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

write.csv(etas, file.path(outDatDir, "3tc_etas.csv"), 
          row.names = FALSE, quote = FALSE)


# QC: 
# Looks like there're a few individuals with very small Ka and Kg, e.g., -3
etas %>% filter(ETA1 < -0.5)
# SUBJECT_ID ID EX_FLAG       ETA1      ETA2          ETA3          ETA4          ETA5       ETA6          ETA7       ETA8
# 1  UMN-034-C 34       0 -0.9168289 0.2585738 -1.238254e-07 -1.754306e-07  0.0129696390  1.5584366 -1.470708e-10 -0.6969728
# 2  UMN-042-C 42       0 -2.7524215 0.1782121 -5.148057e-07 -2.281907e-07 -0.0042844816 -0.2339262 -1.761764e-06 -2.4123175
# 3  UMN-043-C 43       0 -1.5023319 0.3340857  1.942001e-07  1.134508e-07  0.0004811563 -0.1413911  4.635235e-06 -1.4050359

# # Their concentrations are higher than others...what could be the reasons...
obs_final %>% 
   mutate(group = ifelse(ID %in% c(34, 42, 43), 1, 0)) %>% 
   group_by(group, cmt) %>% 
   summarise(meanDV = mean(DV))
# 
 obs_final %>% 
   mutate(group = ifelse(ID %in% c(34, 42, 43), 1, 0)) %>% 
   mutate(group = as.factor(group)) %>% 
   ggplot()+
   geom_point(aes(x=time, y=DV, color=group))+
   facet_wrap(~cmt, ncol=2)+
   scale_x_continuous(limits = c(0, 48))+
   scale_y_log10()


