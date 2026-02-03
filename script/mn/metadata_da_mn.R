
# Library
library(tidyverse)
library(here)
library(readxl)
library(patchwork)
library(lubridate)
library(readr)

# Dir 
sourceDataDir <- here("data/source/mn-biopsy-study")
outDataDir <- here("data/derived/mn")

# Load data ---------------------------------------------------------------

# Function to load all sheets of an excel file
read_data <- function(path, skip){
  sheets_name <- excel_sheets(path)
  data <- map(sheets_name, ~read_xlsx(path, sheet=.x, col_names=TRUE))
  return(data)
}

# Load metadata
metadata <- read_csv(file.path(sourceDataDir, "CTSI_MN_master.csv"))
metadata <- metadata %>%
  rename(
    subject_id = `Patient ID`,
    tfv_p   = `TFV.plasma`,
    ftc_p   = `FTC.plasma`,
    tfvdp_t = `tfvdp.fmol.g`,
    ftctp_t = `ftctp.fmol.g`,
    tfvdp_c = `pbmc.tfvdp`,
    ftctp_c = `pbmc.ftctp`,
    datp_t  = `datp.fmol/g`,
    dctp_t  = `dctp.fmol.g`)
# Missing WT: Subject 002- 87.5 kg, Subject 003- 91.4 kg
metadata <- metadata %>%
  mutate(
    Weight = case_when(
      subject_id == 2 ~ 87.5,
      subject_id == 3 ~ 91.4,
      TRUE ~ Weight)
  )

# Load plasma pk raw data
pconc <- read_xlsx(file.path(sourceDataDir, 
                             "Modified CPAC 1320 Nicol Clinical and Explant Study.xlsx"), 
                   skip=2, sheet="Plasma", col_names=TRUE) 
pconc <- pconc %>%
  mutate(subject_id = as.numeric(`Sample ID`))

# Load tissue pk raw data
tconc <- read_xlsx(file.path(sourceDataDir, 
                             "Modified CPAC 1320 Nicol Clinical and Explant Study.xlsx"), 
                   skip=2, sheet="Tissue", col_names=TRUE)
tconc <- tconc %>%
  mutate(subject_id = as.numeric(`Sample ID`))

# Load PBMC pk raw data
cconc <- read_xlsx(file.path(sourceDataDir, 
                             "Modified CPAC 1320 Nicol Clinical and Explant Study.xlsx"), 
                   skip=2, sheet="Cells", col_names=TRUE)
cconc <- cconc %>%
  mutate(subject_id = as.numeric(`Sample ID`))

# Add sampling time and LLOQs ---------------------------------------------

### PLASMA ---

# Grab sampling time and LLOQ for plasma concentrations (ng/mL)
pconc_pk_time <- pconc %>% 
  mutate(
    dose_datetime = as.POSIXct(
      paste(`Date of Preceding Dose`, format(`Time of Preceding Dose (24hr clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    collect_datetime = as.POSIXct(
      paste(`Date of Sample Collection Post-Dose`, format(`Time of Sample Collection Post-Dose (24hr Clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    pconc_pk_time = as.numeric(difftime(collect_datetime, dose_datetime, units = "hours")),
    plloq = parse_number(`Sample Specific LLOQ`)) %>%  
  dplyr::select(subject_id,
                pconc_pk_time,
                plloq)

# Grab subject ID for plasma concentrations being "BLQ" 
any(grepl("BLQ"           , pconc$`TFV Concentration (ng/mL)`)) # Y
any(grepl("BLQ"           , pconc$`FTC Concentration (ng/mL)`)) # Y

# Pull subject id if he/she has plasma TFV BLQ
pconc_tfv_blq <- pconc %>% 
  filter(`TFV Concentration (ng/mL)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has plasma FTC BLQ
pconc_ftc_blq <- pconc %>% 
  filter(`FTC Concentration (ng/mL)` == "BLQ") %>% 
  pull(subject_id)

### TISSUE ---

# Grab sampling time and LLOQ for tissue concentrations (fmol/g)
tconc_pk_time <- tconc %>% 
  mutate(
    dose_datetime = as.POSIXct(
      paste(`Date of Preceding Dose`, format(`Time of Preceding Dose (24hr clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    collect_datetime = as.POSIXct(
      paste(`Date of Sample Collection Post-Dose`, format(`Time of Sample Collection Post-Dose (24hr Clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    tconc_pk_time = as.numeric(difftime(collect_datetime, dose_datetime, units = "hours"))) %>%  
  dplyr::select(subject_id,
                tconc_pk_time,
                tfvdp_t_lloq  = `TFVdp Sample Specific LLOQ (fmol/g)`, 
                ftctp_t_lloq  = `FTCtp Sample Specific LLOQ (fmol/g)`, 
                datp_t_lloq   = `dATP Sample Specific LLOQ (fmol/g)`, 
                dctp_t_lloq   = `dCTP Sample Specific LLOQ (fmol/g)`)

# Grab subject ID for tissue concentrations being "BLQ" 
any(grepl("BLQ"           , tconc$`TFVdp Concentration (fmol/g)`)) # Y
any(grepl("BLQ"           , tconc$`FTCtp Concentration (fmol/g)`)) # Y
any(grepl("BLQ"           , tconc$`dATP Concentration (fmol/g)`))  # Y
any(grepl("BLQ"           , tconc$`dCTP Concentration (fmol/g)`))  # N

# Pull subject id if he/she has cervical tissue TFVdp BLQ
tconc_tfvdp_blq <- tconc %>% 
  filter(`TFVdp Concentration (fmol/g)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has cervical tissue FTCtp BLQ
tconc_ftcdp_blq <- tconc %>% 
  filter(`FTCtp Concentration (fmol/g)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has cervical tissue dATP BLQ
tconc_datp_blq <- tconc %>% 
  filter(`dATP Concentration (fmol/g)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has cervical tissue dCTP BLQ
tconc_dctp_blq <- tconc %>% 
  filter(`dCTP Concentration (fmol/g)` == "BLQ") %>% 
  pull(subject_id)

### CELL ---

# Grab sampling time and LLOQ for cell concentrations (fmol/millon cells)
cconc_pk_time <- cconc %>% 
  mutate(
    dose_datetime = as.POSIXct(
      paste(`Date of Preceding Dose`, format(`Time of Preceding Dose (24hr clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    collect_datetime = as.POSIXct(
      paste(`Date of Sample Collection Post-Dose`, format(`Time of Sample Collection Post-Dose (24hr Clock)`, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"),
    cconc_pk_time = as.numeric(difftime(collect_datetime, dose_datetime, units = "hours"))) %>%  
  dplyr::select(subject_id, 
                cconc_pk_time,
                # Add missing dATP and dCTP conc in megadata
                datp_c        = `dATP Concentration (fmol/million cells)`,
                dctp_c        = `dCTP Concentration (fmol/million cells)`,
                tfvdp_c_lloq  = `TFVdp Sample Specific LLOQ (fmol/million cells)`, 
                ftctp_c_lloq  = `FTCtp Sample Specific LLOQ (fmol/million cells)`, 
                datp_c_lloq   = `dATP Sample Specific LLOQ (fmol/million cells)`, 
                dctp_c_lloq   = `dCTP Sample Specific LLOQ (fmol/million cells)`)

# Grab subject ID for cell concentrations being "BLQ" or "Not reportable"
any(grepl("BLQ"           , cconc$`TFVdp Concentration (fmol/million cells)`)) # Y
any(grepl("BLQ"           , cconc$`FTCtp Concentration (fmol/million cells)`)) # Y
any(grepl("BLQ"           , cconc$`dATP Concentration (fmol/million cells)`))  # N
any(grepl("BLQ"           , cconc$`dCTP Concentration (fmol/million cells)`))  # N 

# Pull subject id if he/she has PBMC TFVdp BLQ
cconc_tfvdp_blq <- cconc %>% 
  filter(`TFVdp Concentration (fmol/million cells)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has PBMC FTCtp BLQ
cconc_ftcdp_blq <- cconc %>% 
  filter(`FTCtp Concentration (fmol/million cells)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has PBMC dATP BLQ/NR
cconc_datp_blq <- cconc %>% 
  filter(`dATP Concentration (fmol/million cells)` == "BLQ") %>% 
  pull(subject_id)

# Pull subject id if he/she has PBMC dCTP BLQ/NR
cconc_dctp_blq <- cconc %>% 
  filter(`dCTP Concentration (fmol/million cells)` == "BLQ") %>% 
  pull(subject_id)

# Merge sampling time and LLOQ in metadata --------------------------------

# Merge plasma PK sampling time
metadata1 <- metadata %>% 
  left_join(pconc_pk_time) %>%
  # Create a BLQ column for each plasma analyte
  mutate(tfv_p_blq = ifelse(tfv_p == "BLQ" | subject_id %in% pconc_tfv_blq, 1, 0),
         ftc_p_blq = ifelse(ftc_p == "BLQ" | subject_id %in% pconc_ftc_blq, 1, 0)
  ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(tfv_p = ifelse(tfv_p == "BLQ" | subject_id %in% pconc_tfv_blq, plloq/2, as.numeric(tfv_p)), 
         ftc_p = ifelse(ftc_p == "BLQ" | subject_id %in% pconc_ftc_blq, plloq/2, as.numeric(ftc_p))
  ) 

# metadata1 %>% dplyr::select(`subject_id`, contains("_p")) %>% view()

# Merge tissue PK sampling time
metadata2 <- metadata1 %>% 
  left_join(tconc_pk_time) %>%  # dplyr::select(contains("_t")) %>% view()
  # Create a BLQ column for each cervical tissue analyte
  mutate(tfvdp_t_blq  = ifelse(tfvdp_t == "BLQ" | subject_id %in% tconc_tfvdp_blq, 1, 0),
         ftctp_t_blq  = ifelse(ftctp_t == "BLQ" | subject_id %in% tconc_ftcdp_blq, 1, 0), 
         datp_t_blq   = ifelse(datp_t  == "BLQ" | subject_id %in% tconc_datp_blq , 1, 0),
         dctp_t_blq   = ifelse(dctp_t  == "BLQ" | subject_id %in% tconc_dctp_blq , 1, 0) 
  ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(tfvdp_t      = ifelse(tfvdp_t == "BLQ" | subject_id %in% tconc_tfvdp_blq, tfvdp_t_lloq/2, as.numeric(tfvdp_t)), 
         ftctp_t      = ifelse(ftctp_t == "BLQ" | subject_id %in% tconc_ftcdp_blq, ftctp_t_lloq/2, as.numeric(ftctp_t)), 
         datp_t       = ifelse(datp_t  == "BLQ" | subject_id %in% tconc_datp_blq , datp_t_lloq/2 , as.numeric(datp_t)), 
         dctp_t       = ifelse(dctp_t  == "BLQ" | subject_id %in% tconc_dctp_blq , dctp_t_lloq/2 , as.numeric(dctp_t))
  ) 

# metadata2 %>% dplyr::select(`subject_id`, contains("_t")) %>% view()

# Merge cell PK sampling time
metadata3 <- metadata2 %>% 
  left_join(cconc_pk_time) %>% 
  # Create a BLQ column for each PBMC analyte
  mutate(tfvdp_c_blq  = ifelse(tfvdp_c == "BLQ" | subject_id %in% cconc_tfvdp_blq, 1, 0),
         ftctp_c_blq  = ifelse(ftctp_c == "BLQ" | subject_id %in% cconc_ftcdp_blq, 1, 0), 
         datp_c_blq   = ifelse(datp_c  == "BLQ" | subject_id %in% cconc_datp_blq , 1, 0),
         dctp_c_blq   = ifelse(dctp_c  == "BLQ" | subject_id %in% cconc_dctp_blq , 1, 0) 
  ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(tfvdp_c      = ifelse(tfvdp_c == "BLQ" | subject_id %in% cconc_tfvdp_blq, tfvdp_c_lloq/2, as.numeric(tfvdp_c)), 
         ftctp_c      = ifelse(ftctp_c == "BLQ" | subject_id %in% cconc_ftcdp_blq, ftctp_c_lloq/2, as.numeric(ftctp_c)), 
         datp_c       = ifelse(datp_c  == "BLQ" | subject_id %in% cconc_datp_blq , datp_c_lloq/2 , as.numeric(datp_c)), 
         dctp_c       = ifelse(dctp_c  == "BLQ" | subject_id %in% cconc_dctp_blq , dctp_c_lloq/2 , as.numeric(dctp_c))
  ) 

# Quick plot --------------------------------------------------------------

# Plot function
plot_func <- function(.x=pconc_pk_time, .y, .blq,
                      .xlab="Time after dose (hours)", .ylab){
  p <- metadata3 %>% 
    mutate(`BLQ` = ifelse({{.blq}}==1, "Yes", "No")) %>% 
    ggplot()+geom_point(aes(x={{.x}}, y={{.y}}, color=`BLQ`))+
    xlab(.xlab)+ylab(.ylab)+
    theme_bw()#+theme(legend.position = "bottom")
  return(p)
}

# Plasma PK
p1 <- plot_func(.y=`tfv_p` , .blq=`tfv_p_blq` , .ylab="TFV concentration (ng/mL)")
p2 <- plot_func(.y=`ftc_p`, .blq=`ftc_p_blq`, .ylab="FTC concentration (ng/mL)")

p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata3, c("tfv_p//TFV concentration (ng/mL)", 
                                "ftc_p//FTC concentration (ng/mL)"))

# Tissue PK
p3 <- plot_func(.y=`tfvdp_t` , .blq=`tfvdp_t_blq` , .ylab="TFV-dp \nconcentration (fmol/g)")
p4 <- plot_func(.y=`ftctp_t`, .blq=`ftctp_t_blq`, .ylab="FTC-tp \nconcentration (fmol/g)")
p5 <- plot_func(.y=`datp_t`  , .blq=`datp_t_blq`  , .ylab="dATP \nconcentration (fmol/g)"  )
p6 <- plot_func(.y=`dctp_t`  , .blq=`dctp_t_blq`  , .ylab="dCTP \nconcentration (fmol/g)"  )

p3 + p4 + p5 + p6 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata3, c("tfvdp_t//TFV-dp \nconcentration (fmol/g)", 
                                "ftctp_t//FTC-tp \nconcentration (fmol/g)", 
                                "datp_t//dATP \nconcentration (fmol/g)", 
                                "dctp_t//dCTP \nconcentration (fmol/g)"))

# Cell PK
p7  <- plot_func(.y=`tfvdp_c` , .blq=`tfvdp_c_blq` , 
                 .ylab="TFV-dp concentration \n(fmol/million cells)")
p8  <- plot_func(.y=`ftctp_c` , .blq=`ftctp_c_blq` , 
                 .ylab="FTC-tp concentration \n(fmol/million cells)")
p9  <- plot_func(.y=`datp_c`  , .blq=`datp_c_blq`  , 
                 .ylab="dATP concentration \n(fmol/million cells)")
p10 <- plot_func(.y=`dctp_c`  , .blq=`dctp_c_blq`  , 
                 .ylab="dCTP concentration \n(fmol/million cells)")

p7 + p8 + p9 + p10 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata3, c("tfvdp_c//TFV-dp \nconcentration (fmol/million cells)", 
                                "ftctp_c//FTC-tp \nconcentration (fmol/million cells)", 
                                "datp_c//dATP \nconcentration (fmol/million cells)", 
                                "dctp_c//dCTP \nconcentration (fmol/million cells)"))

# Separate metadata -------------------------------------------------------

metadata3 <- metadata3 %>%
  rename_with(~ gsub("\\.", "_", tolower(.x)))
names(metadata3) # 55 variables

# Create vectors of variable names for each subset
covdat_var <- c("subject_id"             ,"age"                    ,"menopasual_status"  ,   
                "hormonal_therapy"       ,"tdf_taf"                ,"group_meno"         ,  
                "group_hor"              ,"ethinicity"             ,"height"             ,   
                "race"                   ,"weight"                 ,"regimen"            ,   
                "evg"                    ,"drv"                    ,"efv"                ,   
                "dtg"                    ,"bic"                    ,"atv"                ,   
                "rtv"                    ,"etr"                    ) # 20 variables

pdat_var <- c("subject_id"  , "pconc_pk_time", 
              "tfv_p"       , "ftc_p"      , 
              "plloq"       , 
              "tfv_p_blq"   , "ftc_p_blq") # 7 variables

tdat_var <- c("subject_id"  , "tconc_pk_time",  
              "tfvdp_t"     , "ftctp_t"     , "datp_t"     , "dctp_t"     , 
              "tfvdp_t_lloq", "ftctp_t_lloq", "datp_t_lloq", "dctp_t_lloq",      
              "tfvdp_t_blq" , "ftctp_t_blq" , "datp_t_blq" , "dctp_t_blq") # 14 variables

cdat_var <- c("subject_id"  , "cconc_pk_time", 
              "tfvdp_c"     , "ftctp_c"      , "datp_c"     , "dctp_c"     , 
              "tfvdp_c_lloq", "ftctp_c_lloq" , "datp_c_lloq", "dctp_c_lloq", 
              "tfvdp_c_blq" , "ftctp_c_blq"  , "datp_c_blq" , "dctp_c_blq") # 14 variables

# Separate `metadata3`
covdat <- metadata3 %>% dplyr::select(all_of(covdat_var))
pdat   <- metadata3 %>% dplyr::select(all_of(pdat_var)  )
tdat   <- metadata3 %>% dplyr::select(all_of(tdat_var)  )
cdat   <- metadata3 %>% dplyr::select(all_of(cdat_var)  )

#' Reshape `pdat`
pdat1 <- pdat %>% dplyr::select(subject_id, pconc_pk_time, plloq, 
                                tfv_p, ftc_p) %>% 
  pivot_longer(cols=c("tfv_p", "ftc_p"), 
               names_to="flag", values_to="dv")

pdat2 <- pdat %>% dplyr::select(subject_id, pconc_pk_time, plloq, 
                                tfv_p_blq, ftc_p_blq) %>% 
  pivot_longer(cols=c("tfv_p_blq", "ftc_p_blq"), 
               names_to="flag", values_to="blq") %>% 
  mutate(flag=gsub("_blq", "", flag))

pdat_final <- left_join(pdat1, pdat2) %>% 
  rename(time=pconc_pk_time, lloq=plloq) %>% 
  relocate(lloq, .after = dv)

#' Reshape `tdat`
tdat1 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t   , ftctp_t, datp_t, dctp_t) %>% 
  pivot_longer(cols=c("tfvdp_t", "ftctp_t", "datp_t", "dctp_t"), 
               names_to="flag", values_to="dv")

tdat2 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t_lloq, ftctp_t_lloq, datp_t_lloq, dctp_t_lloq) %>% 
  pivot_longer(cols=c("tfvdp_t_lloq", "ftctp_t_lloq", "datp_t_lloq", "dctp_t_lloq"), 
               names_to="flag", values_to="lloq") %>% 
  mutate(flag=gsub("_lloq", "", flag))

tdat3 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t_blq, ftctp_t_blq, datp_t_blq, dctp_t_blq) %>% 
  pivot_longer(cols=c("tfvdp_t_blq" , "ftctp_t_blq" , "datp_t_blq" , "dctp_t_blq"), 
               names_to="flag", values_to="blq") %>% 
  mutate(flag=gsub("_blq", "", flag))

tdat_final <- tdat1 %>% 
  left_join(tdat2) %>% 
  left_join(tdat3) %>% 
  rename(time=tconc_pk_time)

#' Reshape `cdat`
cdat1 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c   , ftctp_c, datp_c, dctp_c) %>% 
  pivot_longer(cols=c("tfvdp_c", "ftctp_c", "datp_c", "dctp_c"), 
               names_to="flag", values_to="dv")

cdat2 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c_lloq, ftctp_c_lloq, datp_c_lloq, dctp_c_lloq) %>% 
  pivot_longer(cols=c("tfvdp_c_lloq", "ftctp_c_lloq", "datp_c_lloq", "dctp_c_lloq"), 
               names_to="flag", values_to="lloq") %>% 
  mutate(flag=gsub("_lloq", "", flag))

cdat3 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c_blq, ftctp_c_blq, datp_c_blq, dctp_c_blq) %>% 
  pivot_longer(cols=c("tfvdp_c_blq" , "ftctp_c_blq" , "datp_c_blq" , "dctp_c_blq"), 
               names_to="flag", values_to="blq") %>% 
  mutate(flag=gsub("_blq", "", flag))

cdat_final <- cdat1 %>% 
  left_join(cdat2) %>% 
  left_join(cdat3) %>% 
  rename(time=cconc_pk_time)

# Combine all pk data -----------------------------------------------------

pkdat <- rbind(pdat_final, tdat_final, cdat_final) %>% 
  arrange(subject_id, time)

# NAs in concentrations
map(names(pkdat), function(.x){
  .x <- sym(.x)
  pkdat %>% count(is.na(!!.x))
})

# All NAs are in cell and tissue PK data
pkdat %>% count(is.na(dv), flag)
pkdat %>% filter(is.na(dv)) %>% count(subject_id)
pkdat %>% filter(is.na(dv)) %>% count(subject_id, flag)
# Do not have NA value

# output ------------------------------------------------------------------

write.csv(pkdat , file.path(outDataDir, "pkdat.csv") , row.names=FALSE, quote=FALSE)
write.csv(covdat, file.path(outDataDir, "covdat.csv"), row.names=FALSE, quote=FALSE)

