
# Library and settings ----------------------------------------------------

# Load library
library(here)
library(tidyverse)
library(mrgsolve)
library(patchwork)
library(mrggsave)
library(mrgmisc)
library(pmtables)
library(magrittr)
library(yspec)

# Directories
datDir <- here("data/derived/tfv")
outFigDir <- here("deliv/figure/tfv/ebesim")
if (!dir.exists(outFigDir)) dir.create(outFigDir)
outTabDir <- here("deliv/table/tfv/ebesim")
if (!dir.exists(outTabDir)) dir.create(outTabDir)

# Global options
options(mrgsolve.project=here("model/mrgsolve"), 
        mrggsave.dir = outFigDir,
        pmtables.dir = outTabDir, 
        mrggsave.dev = "pdf", 
        mrg.script = "ebesim.R")

theme_set(theme_bw())

# Load EBEs ---------------------------------------------------------------

etas <- read.csv(file.path(datDir, "tfv_etas.csv"))

# # Simulate 50 doses Q.D. regimen
# etas <- etas %>% 
#   mutate(amt=((300/1000)/635.52)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
#          time=0, cmt=1, evid=1, ii=24, addl=49)

# Assemble doses

# a.	2-1-1 (currently approved in Europe for MSM)
#     (1) double dose 2-24 hours before exposure
#     (2) a dose 24 hours later
#     (3) a dose 48 hours later 

dose1 <- data.frame(ID = unique(etas$ID), 
                    amt=((300*2/1000)/635.52)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
                    time=0, cmt=1, evid=1, ii=0, addl=0) %>% 
  bind_rows(
    data.frame(ID = unique(etas$ID),
               amt=((300/1000)/635.52)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
               time=24, cmt=1, evid=1, ii=24, addl=1)) %>% 
  arrange(ID, time) %>% 
  left_join(etas)

# b.	2-1-1-1-1-1 (total of 7 doses) “expert opinion” for vaginal protection
#     (1) double dose 2-24 hours before exposure
#     (2) a dose 24 hours later
#     (3) a dose 48 hours later 
#     (2) a dose 72 hours later
#     (3) a dose 96 hours later 
#     (2) a dose 120 hours later

dose2 <- data.frame(ID = unique(etas$ID), 
                    amt=((300*2/1000)/635.52)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
                    time=0, cmt=1, evid=1, ii=0, addl=0) %>% 
  bind_rows(
    data.frame(ID = unique(etas$ID), 
               amt=((300/1000)/635.52)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
               time=24, cmt=1, evid=1, ii=24, addl=4)) %>% 
  arrange(ID, time) %>% 
  left_join(etas)

# Load mrgsolve model -----------------------------------------------------

mod <- mread("tdf.cpp")
mod <- update(mod, atol=1e-60, delta=0.1,end=240)
mod <- mod %>% zero_re() # Determinstic simulation, no random effects

mod
param(mod)
revar(mod)

# EBE simulations ---------------------------------------------------------

out1 <- mod %>% 
  data_set(dose1) %>%                # Simulate using EBEs
  mrgsim(obsonly=TRUE,               # Keep observation records only
         etasrc="data.all",          # Search for ETAs in input data
         recover="SUBJECT_ID",       # Recover subject id for merging
         output="df")                # Output simulation results as data frame

out2 <- mod %>% 
  data_set(dose2) %>%                # Simulate using EBEs
  mrgsim(obsonly=TRUE,               # Keep observation records only
         etasrc="data.all",          # Search for ETAs in input data
         recover="SUBJECT_ID",       # Recover subject id for merging
         output="df")                # Output simulation results as data frame

# Load observations -------------------------------------------------------

# Load observations
pkdat <- read.csv(file.path("data/derived", "pkdat.csv"))

pkdat %>% count(flag)
pkdat %>% count(flag, blq)
pkdat %>% count(flag, is.na(dv))

# Subset data
pkdat1 <- pkdat %>% filter(flag %in% c("tfv_p", "tfvdp_t", "datp_t"))

sum <- pkdat1 %>% group_by(flag) %>% 
  summarize(mean = mean(dv, na.rm = TRUE), 
            sd = sd(dv, na.rm = TRUE))

# Process dATP subset
pkdat1 %>% filter(flag=="datp_t"&is.na(dv)) # One dATP as NA 
pkdat_datp <- pkdat1 %>% filter(flag == "datp_t")
pkdat_datp %>% count(blq)
pkdat_datp <- pkdat_datp %>% dplyr::select(subject_id, dv)
pkdat_datp <- pkdat_datp %>% 
  mutate(dv = ifelse(is.na(dv), median(dv, na.rm=TRUE), dv)) %>% 
  rename(datp_t=dv)

# Trim output for plotting
unit_convert <- function(out){
  
  out_plot <- out %>% 
    mutate(tfv_p = (C2/1e15)*287.216*1e6, 
           # fmol/L to ng/mL ( (fmol/L / fmol/mol) * g/mol * 1e6)
           tfvdp_t = C8/1e3) %>% 
    # fmol/L to fmol/mL (fmol/L / mL/L)
    rename(subject_id=SUBJECT_ID) %>% 
    left_join(pkdat_datp) %>% 
    mutate(ratio = tfvdp_t / datp_t)
  
  return(out_plot)
}

out1 <- unit_convert(out1)
out2 <- unit_convert(out2)

# Make a few test plots ---------------------------------------------------

test_plots <- function(out, logy=FALSE){
  
  p1 <- ggplot()+
    geom_line(data=out, aes(x=time, y=tfv_p, group=ID))+
    geom_hline(data=sum[sum$flag == "tfv_p",], 
               aes(yintercept = mean), linetype=2)+
    geom_rect(data=sum[sum$flag == "tfv_p",], 
              aes(xmin=-Inf, xmax=Inf, 
                  ymin=mean-sd, ymax=mean+sd), 
              fill="gray20", alpha=0.25)+
    xlab("Time (hours)")+
    ylab("Plasma TFV concentrations (ng/mL)"); # p1
  
  p2 <- ggplot()+
    geom_line(data=out, aes(x=time, y=tfvdp_t, group=ID))+
    geom_hline(data=sum[sum$flag == "tfvdp_t",], 
               aes(yintercept = mean), linetype=2)+
    geom_rect(data=sum[sum$flag == "tfvdp_t",], 
              aes(xmin=-Inf, xmax=Inf, 
                  ymin=mean-sd, ymax=mean+sd), 
              fill="gray20", alpha=0.25)+
    xlab("Time (hours)")+
    ylab("Cervical tissue TFVdp concentrations (fmol/mL)"); # p2
  
  if(logy){
    p1 <- p1 + scale_y_log10()
    p2 <- p2 + scale_y_log10()}
  
  p <- p1 + p2 + plot_layout(ncol=1)
  
  return(p)
}

test_plots(out1)
test_plots(out2)
test_plots(out1, logy=TRUE)
test_plots(out2, logy=TRUE)

# Plot TFVdp:dATP ratio in tissue -----------------------------------------

ratio_plot <- function(out, logy=FALSE){
  
  p <- ggplot()+
    geom_line(data=out, aes(x=time,y=ratio, group=ID))+
    geom_hline(aes(yintercept = 0.29), linetype=2, color="blue")+
    # From MelanieN, 0.29 in CD4+ T cell 
    geom_hline(aes(yintercept = 0.086), linetype=2, color="red")+
    # From MelanieN, 0.086 in TZM-bl cell 
    xlab("Time (hours)")+
    ylab("Cervical tissue TFVdp:dATP")
  
  if(logy) p <- p + scale_y_log10()
  
  return(p)
} 

p1 <- ratio_plot(out1)
p2 <- ratio_plot(out2)
p1_log <- ratio_plot(out1, logy=TRUE)
p2_log <- ratio_plot(out2, logy=TRUE)

mrggsave(list(p1, p1_log), stem = "ratio_d1", width=6, height=5)
mrggsave(list(p2, p2_log), stem = "ratio_d2", width=6, height=5)

# Plot against TFVdp thresholds -------------------------------------------

thres_plot <- function(out, logy=FALSE){
  
  p <- ggplot()+
    geom_line(data=out, aes(x=time,y=tfvdp_t, group=ID))+
    # geom_hline(aes(yintercept = 716000), linetype=2, color="blue")+
    # # From MelanieN, 716 fmol/mL, EC50 in vaginal explants (Nicol et al. 2015)
    geom_hline(aes(yintercept = 10700), linetype=2, color="red")+
    # From MelanieN, 10.7 fmol/mL, EC50 in cervical explants (Nicol unpublished)
    xlab("Time (hours)")+
    ylab("Cervical tissue TFVdp (fmol/mL)")
  
  if(logy) p <- p + scale_y_log10()
  
  return(p)
}

p3 <- thres_plot(out1)
p4 <- thres_plot(out2)
p3_log <- thres_plot(out1, logy=TRUE)
p4_log <- thres_plot(out2, logy=TRUE)

mrggsave(list(p3, p3_log), stem = "thres_d1", width=6, height=5)
mrggsave(list(p4, p4_log), stem = "thres_d2", width=6, height=5)

# Summarize PK metrics ----------------------------------------------------

sumres <- function(out, name){
  
  xx <- out %>% 
    dplyr::select(ID, time, subject_id, tfv_p, tfvdp_t, datp_t, ratio) %>% 
    mutate(t_above_1 = ifelse(ratio > 0.086, 1, 0)) %>%     # From MelanieN, 0.086 in TZM-bl cell 
    mutate(t_above_2 = ifelse(ratio > 0.29, 1, 0)) %>%      # From MelanieN, 0.29 in CD4+ T cell 
    mutate(t_above_3 = ifelse(tfvdp_t > 10700, 1, 0)) %>%   # From MelanieN, 10700 fmol/mL, EC50 in cervical explants (Nicol unpublished)
    group_by(ID) %>% 
    mutate(cmax_p = max(tfv_p)) %>%                         # Compute plasma Cmax 0-240 hours
    mutate(cmax_t = max(tfvdp_t)) %>%                       # Compute tissue Cmax 0-240 hours
    mutate(auc_p = mrgmisc::auc_partial(
      idv = time, dv = tfv_p, range = c(0, 0+24*10))) %>%   # Compute plasma AUC 0-240 hours
    mutate(auc_t = mrgmisc::auc_partial(
      idv = time, dv = tfvdp_t, range = c(0, 0+24*10))) %>% # Compute tissue AUC 0-240 hours
    mutate(cave_p = auc_p / 240) %>%                        # Compute plasma Cave 0-240 hours
    mutate(cave_t = auc_t / 240) %>%                        # Compute tissue Cave 0-240 hours
    mutate(frac_1 = sum(t_above_1)/n()) %>%                 # Compute fraction of time above thres 1
    mutate(frac_2 = sum(t_above_2)/n()) %>%                 # Compute fraction of time above thres 2
    mutate(frac_3 = sum(t_above_3)/n()) %>%                 # Compute fraction of time above thres 3
    ungroup() %>% 
    distinct(ID, auc_p, auc_t, cmax_p, cmax_t, cave_p, cave_t, 
             frac_1, frac_2, frac_3)
  
  xx <- xx %>% rename_at(vars(auc_p:frac_3), ~paste0(.x, "_", name))
  
  return(xx)
}

summ1 <- sumres(out1, "d1")
summ2 <- sumres(out2, "d2")

pkparam <- out1 %>% distinct(ID, Ka, Vc, CLtt, CLvvtp, CLvetp, TVKg)

summ <- left_join(summ1, summ2) %>% left_join(pkparam)

# print(summ, n=50)
# count(summ, frac_1_d1)
# count(summ, frac_2_d1)
# count(summ, frac_3_d1)

# Summarise fraction time above for 24 hours after last dose --------------

# Last dose time
ldtime1 <- 24*2 # 2-1-1
ldtime2 <- 24*5 # 2-1-1-1-1 

sumres2 <- function(out, ldtime, name){
  
  xx <- out %>% filter(time >= ldtime & time < ldtime+24) %>% # Filter up to 24 hours after last dose
    dplyr::select(ID, time, subject_id, tfv_p, tfvdp_t, datp_t, ratio) %>% 
    mutate(t_above_1 = ifelse(ratio > 0.086, 1, 0)) %>%     # From MelanieN, 0.086 in TZM-bl cell 
    mutate(t_above_2 = ifelse(ratio > 0.29, 1, 0)) %>%      # From MelanieN, 0.29 in CD4+ T cell 
    mutate(t_above_3 = ifelse(tfvdp_t > 10700, 1, 0)) %>%   # From MelanieN, 10700 fmol/mL, EC50 in cervical explants (Nicol unpublished)
    group_by(ID) %>% 
    mutate(frac24_1 = sum(t_above_1)/n()) %>%    # Compute fraction of time above thres 1
    mutate(frac24_2 = sum(t_above_2)/n()) %>%    # Compute fraction of time above thres 2
    mutate(frac24_3 = sum(t_above_3)/n()) %>%    # Compute fraction of time above thres 3
    ungroup() %>% 
    distinct(ID, frac24_1, frac24_2, frac24_3) %>% 
    rename_at(vars(frac24_1:frac24_3), ~paste0(.x, "_", name))
  
  return(xx)
}

summ_ld24_dose1 <- sumres2(out=out1, ldtime=ldtime1, name="d1")
summ_ld24_dose2 <- sumres2(out=out2, ldtime=ldtime2, name="d2")

summ_ld24 <- left_join(summ_ld24_dose1, summ_ld24_dose2)

summ <- left_join(summ, summ_ld24)

# Output tfv exposures ----------------------------------------------------

spec <- yspec::ys_load(file=here(datDir, "tfv_expo.yml")) # Load spec
ys_document(spec, type="working", output_dir=datDir) # Generate a pdf specification file
summ <- summ %>% dplyr::select(names(spec)) # reorder data columns based on spec
ys_check(summ, spec) # Check data set before output
write.csv(summ, file=here(datDir, "tfv_expo.csv"), quote=FALSE, row.names=FALSE)

# Plot PK metrics ---------------------------------------------------------

names <- summ %>% dplyr::select(-ID) %>% names()

p_summ <- map(names, function(.var){
  .var <- as.symbol(.var)
  p <- ggplot(summ)+geom_histogram(aes(x=!!.var))
  return(p)})

p_expo_dose1 <- pmplots::pm_grid(p_summ[1:6]  , ncol=2)
p_frac_dose1 <- pmplots::pm_grid(p_summ[7:12] , ncol=3)
p_expo_dose2 <- pmplots::pm_grid(p_summ[13:18], ncol=2)
p_frac_dose2 <- pmplots::pm_grid(p_summ[19:24], ncol=3)

mrggsave(p_expo_dose1, stem = "expo_d1", width=7, height=8)
mrggsave(p_frac_dose1, stem = "frac_d1", width=8, height=6)
mrggsave(p_expo_dose2, stem = "expo_d2", width=7, height=8)
mrggsave(p_frac_dose2, stem = "frac_d2", width=8, height=6)

p_pairs_dose1 <- summ %>% dplyr::select(contains("d1")) %>% 
  pmplots::pairs_plot(., names[grepl("d1", names)])
mrggsave(p_pairs_dose1, stem = "pairs_d1", width=15, height=15)

p_pairs_dose2 <- summ %>% dplyr::select(contains("d2")) %>% 
  pmplots::pairs_plot(., names[grepl("d2", names)])
mrggsave(p_pairs_dose2, stem = "pairs_d2", width=15, height=15)

# Table PK metrics --------------------------------------------------------

clean_tab <- function(dose){
  
  tab_expo_dose_mean <- summ %>% 
    dplyr::select(contains(dose)) %>% 
    summarise(across(everything(), mean)) %>% 
    pivot_longer(cols=everything(), names_to = "name", values_to = "mean")
  
  tab_expo_dose_sd <- summ %>% 
    dplyr::select(contains(dose)) %>% 
    summarise(across(everything(), sd)) %>% 
    pivot_longer(cols=everything(), names_to = "name", values_to = "sd")
  
  tab_expo_dose <- left_join(tab_expo_dose_mean, tab_expo_dose_sd) %>% 
    mutate(across(mean:sd, pmtables::sig)) %>% 
    mutate(value = paste0(mean, " (", sd, ")")) %>% 
    dplyr::select(-mean, -sd) 
  
  expo <- tab_expo_dose %>% filter(!grepl("frac", name)) %>% 
    mutate(name = gsub(paste0("_", dose), "", name)) %>% 
    separate(col=name, into=c("expo", "tissue"), sep="_") %>% 
    pivot_wider(id_cols=expo, names_from=tissue, values_from=value) %>% 
    rename(`Exposures`=expo, `Plasma`=p, `Cervical`=t) %>% 
    mutate(`Unit1`=c("ng*hour/mL", "ng/mL", "ng/mL"), 
           `Unit2`=c("fmol*hour/mL", "fmol/mL", "fmol/mL")) %>% 
    dplyr::select(`Exposures`,`Plasma`,`Unit1`,`Cervical`,`Unit2`) %>% 
    mutate(Exposures = case_when(Exposures == "auc"~"AUC0-240", 
                                 Exposures == "cmax"~"Cmax", 
                                 Exposures == "cave"~"Cave"))
  
  frac <- tab_expo_dose %>% filter(grepl("frac", name)) %>% 
    mutate(name = gsub(paste0("_", dose), "", name)) %>% 
    mutate(index = gsub("frac", "", name)) %>% 
    mutate(index = gsub("24", "", index)) %>% 
    mutate(name = ifelse(grepl("24", name), "ld24", "all")) %>% 
    pivot_wider(id_cols=name, names_from=index, values_from=value) %>% 
    mutate(name = case_when(name=="all"~"Within 240 hours after first dose", 
                            name=="ld24"~"Within 24 hours after last dose")) %>% 
    rename(`Fraction of time above`=name, 
           `Threshold 1`=`_1`, `Threshold 2`=`_2`, `Threshold 3`=`_3`)
  
  return(list(expo, frac))
}

tab_d1 <- clean_tab("d1")
tab_d2 <- clean_tab("d2")

generate_latex_tab <- function(tab, regimen, name){
  
  tab_expo <- tab[[1]] %>%
    st_new() %>%
    st_rename("Unit"="Unit1", "Unit"="Unit2") %>%
    st_files(output = paste0("tab_expo_",name,".tex")) %>%
    st_notes(paste0("Dosing regimen ", regimen)) %>%
    st_caption(paste0("Summary of exposure metrics using ",regimen," TFV dose")) %>%
    stable() %T>%
    # st2report(ntex = 2) %T>%
    stable_save()
  
  tab_frac <- tab[[2]] %>%
    st_new() %>%
    st_files(output = paste0("tab_frac_", name,".tex")) %>%
    st_notes(paste0("Dosing regimen ", regimen)) %>%
    st_notes("Threshold 1: TFV:dATP ratio of 0.086") %>%
    st_notes("Threshold 2: TFV:dATP ratio of 0.29") %>%
    st_notes("Threshold 3: TFV conc of 10,700 fmol/mL") %>%
    st_caption(paste0("Summary of fraction above threshold using ",regimen," TFV dose")) %>%
    stable() %T>%
    # st2report(ntex = 2) %T>%
    stable_save()
  
  return(list(tab_expo, tab_frac))
}

tab1 <- generate_latex_tab(tab=tab_d1, regimen="2-1-1", name="d1")
tab2 <- generate_latex_tab(tab=tab_d2, regimen="2-1-1-1-1-1", name="d2")

# Make a .pdf table
st2report(
  list(tab1[[1]], tab1[[2]], tab2[[1]], tab2[[2]]), 
  ntex = 2, output_dir = outTabDir,
  stem = "tab_expo",
  show_pdf = interactive()#,
  # dry_run = TRUE,
  #stdout = here(outTabDir, "pdflatex.txt")
)
