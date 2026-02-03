
library(mrgsolve)
library(tidyverse)
library(here)
library(readxl)

options(mrgsolve.project=here("model/mrgsolve"))

theme_set(theme_bw())

# # Load tconc data
# sourceDataDir <- here("data/source")
# tconc <- read_xlsx(file.path(sourceDataDir, 
#                              "Modified CPAC 991 Nicol MUJHU Study (1).xlsx"), 
#                    skip=2, sheet="Tissue", col_names=TRUE)
# tconc <- tconc %>% filter(`Species` == "Human") %>% 
#   mutate(across(`TFVdp Conc (ng/mL)`:`TFVdp Conc (fmol/g)`, as.numeric)) %>% 
#   filter(!is.na(`TFVdp Conc (ng/mL)`))
# 
# tconc %>% ggplot(aes(x=`TFVdp Conc (fmol/g)`, y=`TFVdp Conc (ng/mL)`))+geom_point()
# 
# tconc %>% mutate(coeff = `TFVdp Conc (fmol/g)`/`TFVdp Conc (ng/mL)`) %>% pull(coeff)

mod <- mread("tdf.cpp")
mod <- update(mod, atol=1e-60, delta=0.1)

mod
param(mod)
revar(mod)

#' Molecular weight
#' TDF: 635.52 g/mol
#' Ref: https://pubmed.ncbi.nlm.nih.gov/30150483/
#' TFV: 287.216 g/mol
#' Ref: https://pubmed.ncbi.nlm.nih.gov/30150483/
#' TFVdp: 447.17 g/mol 
#' Ref: https://www.medchemexpress.com/tenofovir-diphosphate.html?srsltid=AfmBOopU1rIbm_rzsOIiFNwPf42v8ZmgwLELWXJouQgyx86f0aOfo253

# Test simulation
mod %>% ev(amt=((300/1000)/635.52)*1e15, # (g / (g/mol))*(fmol/mol)
           time=0, cmt=1, evid=1) %>% 
  mrgsim() %>% plot()

mod %>% ev(amt=((300/1000)/635.52)*1e15, # (g / (g/mol))*(fmol/mol)
           time=0, cmt=1, evid=1, ss=1, ii=24) %>% 
  mrgsim() %>% plot()

mod %>% ev(amt=((300/1000)/635.52)*1e15, # (g / (g/mol))*(fmol/mol)
           time=0, cmt=1, evid=1, ss=1, ii=24) %>% 
  mrgsim(nid=100) %>% plot()

# Simulation
withr::with_seed(
  seed=12513,
  out <- map(1:200, function(i){
    sim <- mod %>% 
      ev(amt=((300/1000)/635.52)*1e15, # (g / (g/mol))*(fmol/mol)
         time=0, cmt=1, evid=1, 
         ss=1, ii=24) %>% # QD dosing to steady-state
      mrgsim(obsonly=TRUE, output="df", 
             nid=100)
    sim <- sim %>% mutate(rep=i)
    return(sim)
   }) %>% bind_rows())

out <- out %>% 
  mutate(tfv_p = (Y2/1e15)*287.216*1e6, # fmol/L to ng/mL ( (fmol/L / fmol/mol) * g/mol * 1e6)
         tfvdp_t = Y8/1e3) # fmol/L to fmol/mL (fmol/L / mL/L)

# Draft plots
out %>% filter(rep==1) %>% ggplot(aes(x=time, y=tfv_p, group=ID))+geom_line()
out %>% filter(rep==1) %>% ggplot(aes(x=time, y=tfvdp_t, group=ID))+geom_line()

# Summary

# Load observations
obs <- read.csv(here("data/derived/pkdat.csv")) %>% 
  mutate(blq = ifelse(blq==1, "BLQ", "Non-BLQ"))

tfv_p_obs <- obs %>% filter(flag == "tfv_p", !is.na(dv))
tfvdp_t_obs <- obs %>% filter(flag == "tfvdp_t", !is.na(dv))

# Summarise simulated PK profiles 
sum <- out %>% 
  group_by(rep, time) %>% 
  summarise(q10_tfv_p   = quantile(tfv_p  , 0.1), 
            q50_tfv_p   = quantile(tfv_p  , 0.5), 
            q90_tfv_p   = quantile(tfv_p  , 0.9), 
            q10_tfvdp_t = quantile(tfvdp_t, 0.1), 
            q50_tfvdp_t = quantile(tfvdp_t, 0.5), 
            q90_tfvdp_t = quantile(tfvdp_t, 0.9)) %>% 
  ungroup() %>% 
  group_by(time) %>% 
  summarise(lo_q10_tfv_p = quantile(q10_tfv_p, 0.025), 
            mi_q10_tfv_p = quantile(q10_tfv_p, 0.5), 
            up_q10_tfv_p = quantile(q10_tfv_p, 0.975), 
            lo_q50_tfv_p = quantile(q50_tfv_p, 0.025), 
            mi_q50_tfv_p = quantile(q50_tfv_p, 0.5), 
            up_q50_tfv_p = quantile(q50_tfv_p, 0.975), 
            lo_q90_tfv_p = quantile(q90_tfv_p, 0.025), 
            mi_q90_tfv_p = quantile(q90_tfv_p, 0.5), 
            up_q90_tfv_p = quantile(q90_tfv_p, 0.975), 
            lo_q10_tfvdp_t = quantile(q10_tfvdp_t, 0.025), 
            mi_q10_tfvdp_t = quantile(q10_tfvdp_t, 0.5), 
            up_q10_tfvdp_t = quantile(q10_tfvdp_t, 0.975), 
            lo_q50_tfvdp_t = quantile(q50_tfvdp_t, 0.025), 
            mi_q50_tfvdp_t = quantile(q50_tfvdp_t, 0.5), 
            up_q50_tfvdp_t = quantile(q50_tfvdp_t, 0.975), 
            lo_q90_tfvdp_t = quantile(q90_tfvdp_t, 0.025), 
            mi_q90_tfvdp_t = quantile(q90_tfvdp_t, 0.5), 
            up_q90_tfvdp_t = quantile(q90_tfvdp_t, 0.975)) %>% 
  ungroup()

sum <- sum %>% filter(time != 0) # Remove time 0 

# Plot plasma simulated versus observed profiles
sum %>% ggplot(aes(x=time))+
  geom_line(aes(y=mi_q10_tfv_p), color="navyblue")+
  geom_ribbon(aes(ymin=lo_q10_tfv_p, ymax=up_q10_tfv_p), alpha=0.2, fill="navyblue")+
  geom_line(aes(y=mi_q50_tfv_p), color="gray20")+
  geom_ribbon(aes(ymin=lo_q50_tfv_p, ymax=up_q50_tfv_p), alpha=0.2, fill="gray20")+
  geom_line(aes(y=mi_q90_tfv_p), color="navyblue")+
  geom_ribbon(aes(ymin=lo_q90_tfv_p, ymax=up_q90_tfv_p), alpha=0.2, fill="navyblue")+
  geom_point(data=tfv_p_obs, aes(x=time, y=dv, shape=blq))+
  geom_hline(data=tfv_p_obs, aes(yintercept=lloq), linetype=2)+
  xlab("Time after dose (hours)")+
  ylab("TFV concentrations (ng/mL)")+
  scale_y_log10()+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")

# Plot cervical simulated versus observed profiles
sum %>% ggplot(aes(x=time))+
  geom_line(aes(y=mi_q10_tfvdp_t), color="navyblue")+
  geom_ribbon(aes(ymin=lo_q10_tfvdp_t, ymax=up_q10_tfvdp_t), alpha=0.2, fill="navyblue")+
  geom_line(aes(y=mi_q50_tfvdp_t), color="gray20")+
  geom_ribbon(aes(ymin=lo_q50_tfvdp_t, ymax=up_q50_tfvdp_t), alpha=0.2, fill="gray20")+
  geom_line(aes(y=mi_q90_tfvdp_t), color="navyblue")+
  geom_ribbon(aes(ymin=lo_q90_tfvdp_t, ymax=up_q90_tfvdp_t), alpha=0.2, fill="navyblue")+
  geom_point(data=tfvdp_t_obs, aes(x=time, y=dv, shape=blq))+
  xlab("Time after dose (hours)")+
  ylab("TFVdp concentrations (fmol/g)")+
  scale_y_log10()+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")


  



