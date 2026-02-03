
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
datDir <- here("data/derived/3tc")
outFigDir <- here("deliv/figure/3tc/ebesim")
if (!dir.exists(outFigDir)) dir.create(outFigDir)
outTabDir <- here("deliv/table/3tc/ebesim")
if (!dir.exists(outTabDir)) dir.create(outTabDir)

# Global options
options(mrgsolve.project=here("model/mrgsolve"), 
        mrggsave.dir = outFigDir,
        pmtables.dir = outTabDir, 
        mrggsave.dev = "pdf", 
        mrg.script = "ebesim_3tc.R")

theme_set(theme_bw())

# Load EBEs ---------------------------------------------------------------

etas <- read.csv(file.path(datDir, "3tc_etas.csv"))

# # Simulate 50 doses Q.D. regimen
# etas <- etas %>% 
#   mutate(amt=((300/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
#          time=0, cmt=1, evid=1, ii=24, addl=49)

# Assemble doses

# a.	2-1-1 (currently approved in Europe for MSM)
#     (1) double dose 2-24 hours before exposure
#     (2) a dose 24 hours later
#     (3) a dose 48 hours later 

dose1 <- data.frame(ID = unique(etas$ID), 
                    amt=((300*2/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
                    time=0, cmt=1, evid=1, ii=0, addl=0) %>% 
  bind_rows(
    data.frame(ID = unique(etas$ID),
               amt=((300/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
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
                    amt=((300*2/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
                    time=0, cmt=1, evid=1, ii=0, addl=0) %>% 
  bind_rows(
    data.frame(ID = unique(etas$ID), 
               amt=((300/1000)/229.26)*1e15, # fmol dosing (g / (g/mol))*(fmol/mol)
               time=24, cmt=1, evid=1, ii=24, addl=4)) %>% 
  arrange(ID, time) %>% 
  left_join(etas)

# Load mrgsolve model -----------------------------------------------------

mod <- mread("3tc.cpp")
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
pkdat1 <- pkdat %>% filter(flag %in% c("_3tc_p", "_3tctp_t"))

sum <- pkdat1 %>% group_by(flag) %>% 
  summarize(mean = mean(dv, na.rm = TRUE), 
            sd = sd(dv, na.rm = TRUE))

# Trim output for plotting
unit_convert <- function(out){
  
  out_plot <- out %>% 
    mutate(`_3tc_p` = (C2/1e15)*229.26*1e6, 
           # fmol/L to ng/mL ( (fmol/L / fmol/mol) * g/mol * 1e6)
           `_3tctp_t` = C8/1e3) %>% 
    # fmol/L to fmol/mL (fmol/L / mL/L)
    rename(subject_id=SUBJECT_ID)
  
  return(out_plot)
}

out1 <- unit_convert(out1)
out2 <- unit_convert(out2)

# Make a few test plots ---------------------------------------------------

test_plots <- function(out, logy=FALSE){
  
  p1 <- ggplot()+
    geom_line(data=out, aes(x=time, y=`_3tc_p`, group=ID))+
    geom_hline(data=sum[sum$flag == "_3tc_p",], 
               aes(yintercept = mean), linetype=2)+
    geom_rect(data=sum[sum$flag == "_3tc_p",], 
              aes(xmin=-Inf, xmax=Inf, 
                  ymin=mean-sd, ymax=mean+sd), 
              fill="gray20", alpha=0.25)+
    xlab("Time (hours)")+
    ylab("Plasma 3TC concentrations (ng/mL)"); # p1
  
  p2 <- ggplot()+
    geom_line(data=out, aes(x=time, y=`_3tctp_t`, group=ID))+
    geom_hline(data=sum[sum$flag == "_3tctp_t",], 
               aes(yintercept = mean), linetype=2)+
    geom_rect(data=sum[sum$flag == "_3tctp_t",], 
              aes(xmin=-Inf, xmax=Inf, 
                  ymin=mean-sd, ymax=mean+sd), 
              fill="gray20", alpha=0.25)+
    xlab("Time (hours)")+
    ylab("Cervical tissue 3TCtp concentrations (fmol/mL)"); # p2
  
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

# Summarize PK metrics ----------------------------------------------------

sumres <- function(out, name){
  
  xx <- out %>% 
    dplyr::select(ID, time, subject_id, `_3tc_p`, `_3tctp_t`) %>% 
    group_by(ID) %>% 
    mutate(cmax_p = max(`_3tc_p`)) %>%                         # Compute plasma Cmax 0-240 hours
    mutate(cmax_t = max(`_3tctp_t`)) %>%                       # Compute tissue Cmax 0-240 hours
    mutate(auc_p = mrgmisc::auc_partial(
      idv = time, dv = `_3tc_p`, range = c(0, 0+24*10))) %>%   # Compute plasma AUC 0-240 hours
    mutate(auc_t = mrgmisc::auc_partial(
      idv = time, dv = `_3tctp_t`, range = c(0, 0+24*10))) %>% # Compute tissue AUC 0-240 hours
    mutate(cave_p = auc_p / 240) %>%                        # Compute plasma Cave 0-240 hours
    mutate(cave_t = auc_t / 240) %>%                        # Compute tissue Cave 0-240 hours
    ungroup() %>% 
    distinct(ID, auc_p, auc_t, cmax_p, cmax_t, cave_p, cave_t)
  
  xx <- xx %>% rename_at(vars(auc_p:cave_t), ~paste0(.x, "_", name))
  
  return(xx)
}

summ1 <- sumres(out1, "d1")
summ2 <- sumres(out2, "d2")

pkparam <- out1 %>% distinct(ID, Ka, Vc, CLttvvtp, CLvvtp, CLttve, CLvetp, CLttvr, Kg)

summ <- left_join(summ1, summ2) %>% left_join(pkparam)

# print(summ, n=50)

# Output 3TC exposures ----------------------------------------------------

spec <- yspec::ys_load(file=here(datDir, "3tc_expo.yml")) # Load spec
ys_document(spec, type="working", output_dir=datDir) # Generate a pdf specification file
summ <- summ %>% dplyr::select(names(spec)) # reorder data columns based on spec
ys_check(summ, spec) # Check data set before output
write.csv(summ, file=here(datDir, "3tc_expo.csv"), quote=FALSE, row.names=FALSE)

# Plot PK metrics ---------------------------------------------------------

names <- summ %>% dplyr::select(-ID) %>% names()

p_summ <- map(names, function(.var){
  .var <- as.symbol(.var)
  p <- ggplot(summ)+geom_histogram(aes(x=!!.var))
  return(p)})

p_param <- pmplots::pm_grid(p_summ[1:8]  , ncol=2)
p_expo_dose1 <- pmplots::pm_grid(p_summ[9:14]  , ncol=2)
p_expo_dose2 <- pmplots::pm_grid(p_summ[15:20], ncol=2)

mrggsave(p_param, stem = "parameters", width=7, height=8)
mrggsave(p_expo_dose1, stem = "expo_d1", width=7, height=8)
mrggsave(p_expo_dose2, stem = "expo_d2", width=7, height=8)

p_pairs_dose1 <- summ %>% dplyr::select(contains("d1")) %>% 
  pmplots::pairs_plot(., names[grepl("d1", names)])
mrggsave(p_pairs_dose1, stem = "pairs_d1", width=10, height=10)

p_pairs_dose2 <- summ %>% dplyr::select(contains("d2")) %>% 
  pmplots::pairs_plot(., names[grepl("d2", names)])
mrggsave(p_pairs_dose2, stem = "pairs_d2", width=10, height=10)

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
  return(expo)
}

tab_d1 <- clean_tab("d1")
tab_d2 <- clean_tab("d2")

generate_latex_tab <- function(tab, regimen, name){
  
  tab_expo <- tab %>%
    st_new() %>%
    st_rename("Unit"="Unit1", "Unit"="Unit2") %>%
    st_files(output = paste0("tab_expo_",name,".tex")) %>%
    st_notes(paste0("Dosing regimen ", regimen)) %>%
    st_caption(paste0("Summary of exposure metrics using ",regimen," 3TC dose")) %>%
    stable() %T>%
    # st2report(ntex = 2) %T>%
    stable_save()
  
  return(tab_expo)
}

tab1 <- generate_latex_tab(tab=tab_d1, regimen="2-1-1", name="d1")
tab2 <- generate_latex_tab(tab=tab_d2, regimen="2-1-1-1-1-1", name="d2")

# Make a .pdf table
st2report(
  list(tab1, tab2), 
  ntex = 2, output_dir = outTabDir,
  stem = "tab_expo",
  show_pdf = interactive()#,
  # dry_run = TRUE,
  #stdout = here(outTabDir, "pdflatex.txt")
)
