
# Library
library(tidyverse)
library(here)
library(readxl)
library(patchwork)

# Dir 
sourceDataDir <- here("data/source")
outDataDir <- here("data/derived")

# Load data ---------------------------------------------------------------

# Function to load all sheets of an excel file
read_data <- function(path, skip){
  sheets_name <- excel_sheets(path)
  data <- map(sheets_name, ~read_xlsx(path, sheet=.x, col_names=TRUE))
  return(data)
}

# Load metadata
metadata <- read_xlsx(file.path(sourceDataDir, "metadata 7Feb2019.xlsx"), 
                      sheet="data", col_names=TRUE)

# Load plasma pk raw data
pconc <- read_xlsx(file.path(sourceDataDir, 
                             "Modified CPAC 991 Nicol MUJHU Study_Plasma.xlsx"), 
                   skip=2, sheet="Plasma", col_names=TRUE) %>% 
  # Two typos in `Treatment group (Optional)` found, corrections made
  mutate(`Treatment group (Optional)` = ifelse(`Treatment group (Optional)` == "UMN-014A", 
                                               "UMN-014-A", `Treatment group (Optional)`)) %>% 
  mutate(`Treatment group (Optional)` = ifelse(`Treatment group (Optional)` == "UMN-021B", 
                                               "UMN-021-B", `Treatment group (Optional)`))

# Load tissue pk raw data
tconc <- read_xlsx(file.path(sourceDataDir, 
                             "Modified CPAC 991 Nicol MUJHU Study (1).xlsx"), 
                   skip=2, sheet="Tissue", col_names=TRUE)

# Load PBMC pk raw data
cconc <- read_xlsx(file.path(sourceDataDir, 
                             "PBMC Concentrations from MUJHU study.xlsx"), 
                   skip=2, sheet="Cells", col_names=TRUE)

# Load pk sample time raw data
tdose <- read_csv(file.path(sourceDataDir, 
                            "ModifiersOfTenofovir_time from last dose.csv"))

# Some mismatch in data input detected ------------------------------------

# A few inconsistent NAs apears in `metadata` and tissue PK raw data
# For example, UMN-027-B and UMN-033-B
test1 <- metadata %>% 
  dplyr::select(subject_id, 
                `tfvdp_t_metadata`  = `tfvdp_t` , 
                `_3tctp_t_metadata` = `_3tctp_t`, 
                `datp_t_metadata`   = `datp_t`  , 
                `dctp_t_metadata`   = `dctp_t`)

test2 <- tconc %>% 
  dplyr::select(`subject_id`      = `Treatment group (Optional)`, 
                `tfvdp_t_tconc`   = `TFVdp Conc (fmol/g)`, 
                `_3tctp_t_tconc`  = `3TCtp Conc (fmol/g)`, 
                `datp_t_tconc`    = `dATP Conc (fmol/g)` , 
                `dctp_t_tconc`    = `dCTP Conc (fmol/g)` )

left_join(test1, test2) %>% 
  mutate(across(tfvdp_t_metadata:dctp_t_tconc, as.numeric)) %>% 
  mutate(across(`tfvdp_t_tconc`:`dctp_t_tconc`, ~round(.x, digits=0))) %>% 
  view()

# A few inconsistent NAs apears in `metadata` and PBMC PK raw data
# For example, UMN-027-B, UMN-038-C
test3 <- metadata %>% 
  dplyr::select(subject_id, 
                `tfvdp_c_metadata`  = `tfvdp_c`, 
                `3tctp_c_metadata`  = `3tctp_c`, 
                `datp_c_metadata`   = `datp_c` , 
                `dctp_c_metadata`   = `dctp_c`)

test4 <- cconc %>% 
  dplyr::select(`subject_id`      = `Treatment group (Optional)`, 
                `tfvdp_c_tconc`   = `TFVdp Conc (fmol/million cells)`, 
                `3tctp_c_tconc`   = `3TCtp Conc (fmol/million cells)`, 
                `datp_c_tconc`    = `dATP Conc (fmol/million cells)` , 
                `dctp_c_tconc`    = `dCTP Conc (fmol/million cells)` )

left_join(test3, test4) %>% 
  mutate(across(tfvdp_c_metadata:dctp_c_tconc, as.numeric)) %>% 
  mutate(across(`tfvdp_c_tconc`:`dctp_c_tconc`, ~round(.x, digits=0))) %>% 
  view()

# Add sampling time and LLOQs ---------------------------------------------

### PLASMA ---

# Grab sampling time and LLOQ for plasma concentrations (ng/mL)
pconc_pk_time <- pconc %>% 
  dplyr::select(`subject_id`    = `Treatment group (Optional)`, 
                `pconc_pk_time` = `PK Timepoint (nominal time post-dose)`, 
                `plloq`         = `Sample Specific LLOQ (ng/mL)`)

# Grab subject ID for plasma concentrations being "BLQ" or "Not reportable"
any(grepl("BLQ"           , pconc$`TFV Concentration (ng/mL)`)) # Y
any(grepl("Not reportable", pconc$`TFV Concentration (ng/mL)`)) # N
any(grepl("BLQ"           , pconc$`3TC Concentration (ng/mL)`)) # Y
any(grepl("Not reportable", pconc$`3TC Concentration (ng/mL)`)) # N

# Pull subject id if he/she has plasma TFV BLQ
pconc_tfv_blq <- pconc %>% 
  filter(`TFV Concentration (ng/mL)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has plasma 3TC BLQ
pconc_3tc_blq <- pconc %>% 
  filter(`3TC Concentration (ng/mL)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

### TISSUE ---

# Grab sampling time and LLOQ for tissue concentrations (fmol/g)
tconc_pk_time <- tconc %>% 
  dplyr::select(`subject_id`    = `Treatment group (Optional)`, 
                `tconc_pk_time` = `PK Timepoint (nominal time post-dose)(hours)`, 
                `tfvdp_t_lloq`  = `Sample Specific LLOQ...20`, 
                `_3tctp_t_lloq` = `Sample Specific LLOQ...32`, 
                `datp_t_lloq`   = `Sample Specific LLOQ...26`, 
                `dctp_t_lloq`   = `Sample Specific LLOQ...29`)

# Grab subject ID for tissue concentrations being "BLQ" or "Not reportable"
any(grepl("BLQ"           , tconc$`TFVdp Conc (fmol/g)`)) # Y
any(grepl("Not reportable", tconc$`TFVdp Conc (fmol/g)`)) # N
any(grepl("BLQ"           , tconc$`3TCtp Conc (fmol/g)`)) # Y
any(grepl("Not reportable", tconc$`3TCtp Conc (fmol/g)`)) # N
any(grepl("BLQ"           , tconc$`dATP Conc (fmol/g)`))  # Y
any(grepl("Not reportable", tconc$`dATP Conc (fmol/g)`))  # N
any(grepl("BLQ"           , tconc$`dCTP Conc (fmol/g)`))  # Y
any(grepl("Not reportable", tconc$`dCTP Conc (fmol/g)`))  # N

# Pull subject id if he/she has cervical tissue TFVdp BLQ
tconc_tfvdp_blq <- tconc %>% 
  filter(`TFVdp Conc (fmol/g)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has cervical tissue 3TCtp BLQ
tconc_3tcdp_blq <- tconc %>% 
  filter(`3TCtp Conc (fmol/g)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has cervical tissue dATP BLQ
tconc_datp_blq <- tconc %>% 
  filter(`dATP Conc (fmol/g)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has cervical tissue dCTP BLQ
tconc_dctp_blq <- tconc %>% 
  filter(`dCTP Conc (fmol/g)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)

### CELL ---

# Grab sampling time and LLOQ for cell concentrations (fmol/millon cells)
cconc_pk_time <- cconc %>% 
  dplyr::select(`subject_id`    =`Treatment group (Optional)`, 
                `cconc_pk_time` = `PK Timepoint (nominal time post-dose)(hours)`, 
                `tfvdp_c_lloq`  = `Sample Specific LLOQ (TFVdp fmol/million cells)`, 
                `_3tctp_c_lloq` = `Sample Specific LLOQ (3TCtp fmol/million cells)`, 
                `datp_c_lloq`   = `Sample Specific LLOQ (dATP fmol/million cells)`, 
                `dctp_c_lloq`   = `Sample Specific LLOQ (dCTP fmol/million cells)`)

# Grab subject ID for cell concentrations being "BLQ" or "Not reportable"
any(grepl("BLQ"           , cconc$`TFVdp Conc (fmol/million cells)`)) # Y
any(grepl("Not reportable", cconc$`TFVdp Conc (fmol/million cells)`)) # Y
any(grepl("BLQ"           , cconc$`3TCtp Conc (fmol/million cells)`)) # Y
any(grepl("Not reportable", cconc$`3TCtp Conc (fmol/million cells)`)) # Y
any(grepl("BLQ"           , cconc$`dATP Conc (fmol/million cells)`))  # N
any(grepl("Not reportable", cconc$`dATP Conc (fmol/million cells)`))  # Y
any(grepl("BLQ"           , cconc$`dCTP Conc (fmol/million cells)`))  # N 
any(grepl("Not reportable", cconc$`dCTP Conc (fmol/million cells)`))  # Y

# Pull subject id if he/she has PBMC TFVdp BLQ/NR
cconc_tfvdp_blq <- cconc %>% 
  filter(`TFVdp Conc (fmol/million cells)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)
cconc_tfvdp_nr  <- cconc %>% 
  filter(`TFVdp Conc (fmol/million cells)` == "Not reportable") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has PBMC 3TCtp BLQ/NR
cconc_3tcdp_blq <- cconc %>% 
  filter(`3TCtp Conc (fmol/million cells)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)
cconc_3tcdp_nr  <- cconc %>% 
  filter(`3TCtp Conc (fmol/million cells)` == "Not reportable") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has PBMC dATP BLQ/NR
cconc_datp_blq <- cconc %>% 
  filter(`dATP Conc (fmol/million cells)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)
cconc_datp_nr <- cconc %>% 
  filter(`dATP Conc (fmol/million cells)` == "Not reportable") %>% 
  pull(`Treatment group (Optional)`)

# Pull subject id if he/she has PBMC dCTP BLQ/NR
cconc_dctp_blq <- cconc %>% 
  filter(`dCTP Conc (fmol/million cells)` == "BLQ") %>% 
  pull(`Treatment group (Optional)`)
cconc_dctp_nr <- cconc %>% 
  filter(`dCTP Conc (fmol/million cells)` == "Not reportable") %>% 
  pull(`Treatment group (Optional)`)

# Merge sampling time and LLOQ in metadata --------------------------------

# Merge plasma PK sampling time
metadata1 <- metadata %>% 
  left_join(pconc_pk_time) %>%
  # Create a BLQ column for each plasma analyte
  mutate(`tfv_p_blq`  = ifelse(is.na(`tfv_p`)  & subject_id %in% pconc_tfv_blq, 1, 0),
         `_3tc_p_blq` = ifelse(is.na(`_3tc_p`) & subject_id %in% pconc_3tc_blq, 1, 0)
  ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(`tfv_p`  = ifelse(is.na(`tfv_p`)  & subject_id %in% pconc_tfv_blq, plloq/2, `tfv_p` ), 
         `_3tc_p` = ifelse(is.na(`_3tc_p`) & subject_id %in% pconc_3tc_blq, plloq/2, `_3tc_p`)
         ) 

# metadata1 %>% dplyr::select(`subject_id`, contains("_p")) %>% view()

# Merge tissue PK sampling time
metadata2 <- metadata1 %>% 
  left_join(tconc_pk_time) %>% # dplyr::select(contains("_t")) %>% view()
  mutate(across(tfvdp_t_lloq:dctp_t_lloq, as.numeric)) %>% # view()
  # Create a BLQ column for each cervical tissue analyte
  mutate(`tfvdp_t_blq`  = ifelse(is.na(`tfvdp_t`)  & subject_id %in% tconc_tfvdp_blq, 1, 0),
         `_3tctp_t_blq` = ifelse(is.na(`_3tctp_t`) & subject_id %in% tconc_3tcdp_blq, 1, 0), 
         `datp_t_blq`   = ifelse(is.na(`datp_t`)   & subject_id %in% tconc_datp_blq , 1, 0),
         `dctp_t_blq`   = ifelse(is.na(`dctp_t`)   & subject_id %in% tconc_dctp_blq , 1, 0) 
         ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(`tfvdp_t`      = ifelse(is.na(`tfvdp_t`)  & subject_id %in% tconc_tfvdp_blq, `tfvdp_t_lloq`/2 , `tfvdp_t` ), 
         `_3tctp_t`     = ifelse(is.na(`_3tctp_t`) & subject_id %in% tconc_3tcdp_blq, `_3tctp_t_lloq`/2, `_3tctp_t`), 
         `datp_t`       = ifelse(is.na(`datp_t`)   & subject_id %in% tconc_datp_blq , `datp_t_lloq`/2. , `datp_t`  ), 
         `dctp_t`       = ifelse(is.na(`dctp_t`)   & subject_id %in% tconc_dctp_blq , `dctp_t_lloq`/2. , `dctp_t`  )
         ) 

# metadata2 %>% dplyr::select(`subject_id`, contains("_t")) %>% view()

# Merge cell PK sampling time
metadata3 <- metadata2 %>% 
  left_join(cconc_pk_time) %>% 
  # Create a BLQ column for each PBMC analyte
  mutate(`tfvdp_c_blq`  = ifelse(is.na(`tfvdp_c`)  & subject_id %in% cconc_tfvdp_blq, 1, 0),
         `3tctp_c_blq`  = ifelse(is.na(`3tctp_c`)  & subject_id %in% cconc_3tcdp_blq, 1, 0), 
         `datp_c_blq`   = ifelse(is.na(`datp_c`)   & subject_id %in% cconc_datp_blq , 1, 0), 
         `dctp_c_blq`   = ifelse(is.na(`dctp_c`)   & subject_id %in% cconc_dctp_blq , 1, 0)
  ) %>% 
  # Fill in DV=LLOQ/2 if BLQ
  mutate(`tfvdp_c`      = ifelse(is.na(`tfvdp_c`)  & subject_id %in% cconc_tfvdp_blq, `tfvdp_c_lloq`/2 , `tfvdp_c`), 
         `3tctp_c`      = ifelse(is.na(`3tctp_c`)  & subject_id %in% cconc_3tcdp_blq, `_3tctp_c_lloq`/2, `3tctp_c`), 
         `datp_c`       = ifelse(is.na(`datp_c`)   & subject_id %in% cconc_datp_blq , `datp_c_lloq`/2  , `datp_c` ), 
         `dctp_c`       = ifelse(is.na(`dctp_c`)   & subject_id %in% cconc_dctp_blq , `dctp_c_lloq`/2  , `dctp_c` )
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
p2 <- plot_func(.y=`_3tc_p`, .blq=`_3tc_p_blq`, .ylab="3TC concentration (ng/mL)")

p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata, c("tfv_p//TFV concentration (ng/mL)", 
                                "_3tc_p//3TC concentration (ng/mL)"))

# Tissue PK
p3 <- plot_func(.y=`tfvdp_t` , .blq=`tfvdp_t_blq` , .ylab="TFV-dp \nconcentration (fmol/g)")
p4 <- plot_func(.y=`_3tctp_t`, .blq=`_3tctp_t_blq`, .ylab="3TC-tp \nconcentration (fmol/g)")
p5 <- plot_func(.y=`datp_t`  , .blq=`datp_t_blq`  , .ylab="dATP \nconcentration (fmol/g)"  )
p6 <- plot_func(.y=`dctp_t`  , .blq=`dctp_t_blq`  , .ylab="dCTP \nconcentration (fmol/g)"  )

p3 + p4 + p5 + p6 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata, c("tfvdp_t//TFV-dp \nconcentration (fmol/g)", 
                                "_3tctp_t//3TC-tp \nconcentration (fmol/g)", 
                                "datp_t//dATP \nconcentration (fmol/g)", 
                                "dctp_t//dCTP \nconcentration (fmol/g)"))

# Cell PK
p7  <- plot_func(.y=`tfvdp_c` , .blq=`tfvdp_c_blq` , 
                 .ylab="TFV-dp concentration \n(fmol/million cells)")
p8  <- plot_func(.y=`3tctp_c` , .blq=`3tctp_c_blq` , 
                 .ylab="3TC-tp concentration \n(fmol/million cells)")
p9  <- plot_func(.y=`datp_c`  , .blq=`datp_c_blq`  , 
                 .ylab="dATP concentration \n(fmol/million cells)")
p10 <- plot_func(.y=`dctp_c`  , .blq=`dctp_c_blq`  , 
                 .ylab="dCTP concentration \n(fmol/million cells)")

p7 + p8 + p9 + p10 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pmplots::pairs_plot(metadata, c("tfvdp_c//TFV-dp \nconcentration (fmol/million cells)", 
                                "3tctp_c//3TC-tp \nconcentration (fmol/million cells)", 
                                "datp_c//dATP \nconcentration (fmol/million cells)", 
                                "dctp_c//dCTP \nconcentration (fmol/million cells)"))

# Separate metadata -------------------------------------------------------

names(metadata3) # 49 variables

# Create vectors of variable names for each subset
covdat_var <- c("subject_id"       , "age"              , "marital_status"  , 
                "sexually_active"  , "ammenorheic"      , "LMP_last"        , 
                "weight"           , "height"           , "BMI"             , 
                "urinary_gonorrhea", "urinary_chlamydia", "syphilis_results", 
                "arv_1"            , "E2"               , "P4"              , 
                "contra"           , "hormone") # 17 variables

pdat_var <- c("subject_id"  , "pconc_pk_time", 
               "tfv_p"       , "_3tc_p"      , 
               "plloq"       , 
               "tfv_p_blq"   , "_3tc_p_blq") # 7 variables

tdat_var <- c("subject_id"  , "tconc_pk_time",  
              "tfvdp_t"     , "_3tctp_t"     , "datp_t"     , "dctp_t"     , 
              "tfvdp_t_lloq", "_3tctp_t_lloq", "datp_t_lloq", "dctp_t_lloq",      
              "tfvdp_t_blq" , "_3tctp_t_blq" , "datp_t_blq" , "dctp_t_blq") # 14 variables

cdat_var <- c("subject_id"  , "cconc_pk_time", 
              "tfvdp_c"     , "3tctp_c"      , "datp_c"     , "dctp_c"     , 
              "tfvdp_c_lloq", "_3tctp_c_lloq", "datp_c_lloq", "dctp_c_lloq", 
              "tfvdp_c_blq" , "3tctp_c_blq"  , "datp_c_blq" , "dctp_c_blq") # 14 variables

# Separate `metadata3`
covdat <- metadata3 %>% dplyr::select(all_of(covdat_var))
pdat   <- metadata3 %>% dplyr::select(all_of(pdat_var)  )
tdat   <- metadata3 %>% dplyr::select(all_of(tdat_var)  )
cdat   <- metadata3 %>% dplyr::select(all_of(cdat_var)  )

#' Reshape `pdat`
pdat1 <- pdat %>% dplyr::select(subject_id, pconc_pk_time, plloq, 
                                tfv_p, `_3tc_p`) %>% 
  pivot_longer(cols=c("tfv_p", "_3tc_p"), 
               names_to="flag", values_to="dv")

pdat2 <- pdat %>% dplyr::select(subject_id, pconc_pk_time, plloq, 
                                tfv_p_blq, `_3tc_p_blq`) %>% 
  pivot_longer(cols=c("tfv_p_blq", "_3tc_p_blq"), 
               names_to="flag", values_to="blq") %>% 
  mutate(flag=gsub("_blq", "", flag))

pdat_final <- left_join(pdat1, pdat2) %>% 
  rename(time=pconc_pk_time, lloq=plloq) %>% 
  relocate(lloq, .after = dv)

#' Reshape `tdat`
tdat1 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t   , `_3tctp_t`, datp_t, dctp_t) %>% 
  pivot_longer(cols=c("tfvdp_t", "_3tctp_t", "datp_t", "dctp_t"), 
               names_to="flag", values_to="dv")

tdat2 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t_lloq, `_3tctp_t_lloq`, datp_t_lloq, dctp_t_lloq) %>% 
  pivot_longer(cols=c("tfvdp_t_lloq", "_3tctp_t_lloq", "datp_t_lloq", "dctp_t_lloq"), 
               names_to="flag", values_to="lloq") %>% 
  mutate(flag=gsub("_lloq", "", flag))

tdat3 <- tdat %>% dplyr::select(subject_id, tconc_pk_time, 
                                tfvdp_t_blq, `_3tctp_t_blq`, datp_t_blq, dctp_t_blq) %>% 
  pivot_longer(cols=c("tfvdp_t_blq" , "_3tctp_t_blq" , "datp_t_blq" , "dctp_t_blq"), 
               names_to="flag", values_to="blq") %>% 
  mutate(flag=gsub("_blq", "", flag))

tdat_final <- tdat1 %>% 
  left_join(tdat2) %>% 
  left_join(tdat3) %>% 
  rename(time=tconc_pk_time)

#' Reshape `cdat`
cdat1 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c   , `_3tctp_c`=`3tctp_c`, datp_c, dctp_c) %>% 
  pivot_longer(cols=c("tfvdp_c", "_3tctp_c", "datp_c", "dctp_c"), 
               names_to="flag", values_to="dv")

cdat2 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c_lloq, `_3tctp_c_lloq`, datp_c_lloq, dctp_c_lloq) %>% 
  pivot_longer(cols=c("tfvdp_c_lloq", "_3tctp_c_lloq", "datp_c_lloq", "dctp_c_lloq"), 
               names_to="flag", values_to="lloq") %>% 
  mutate(flag=gsub("_lloq", "", flag))

cdat3 <- cdat %>% 
  dplyr::select(subject_id, cconc_pk_time, 
                tfvdp_c_blq, `_3tctp_c_blq`=`3tctp_c_blq`, datp_c_blq, dctp_c_blq) %>% 
  pivot_longer(cols=c("tfvdp_c_blq" , "_3tctp_c_blq" , "datp_c_blq" , "dctp_c_blq"), 
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
# # A tibble: 6 × 2
# subject_id     n
# <chr>      <int>
# 1 UMN-019-B      4
# 2 UMN-020-C      2
# 3 UMN-021-B      3
# 4 UMN-027-B      2
# 5 UMN-033-B      3
# 6 UMN-038-C      1
pkdat %>% filter(is.na(dv)) %>% count(subject_id, flag)
# # A tibble: 15 × 3
# subject_id flag         n
# <chr>      <chr>    <int>
# 1 UMN-019-B  _3tctp_c     1
# 2 UMN-019-B  datp_c       1
# 3 UMN-019-B  dctp_c       1
# 4 UMN-019-B  tfvdp_c      1
# 5 UMN-020-C  datp_c       1
# 6 UMN-020-C  dctp_c       1
# 7 UMN-021-B  datp_c       1
# 8 UMN-021-B  dctp_c       1
# 9 UMN-021-B  tfvdp_c      1
# 10 UMN-027-B  _3tctp_c     1
# 11 UMN-027-B  tfvdp_t      1
# 12 UMN-033-B  _3tctp_t     1
# 13 UMN-033-B  datp_t       1
# 14 UMN-033-B  dctp_t       1
# 15 UMN-038-C  _3tctp_c     1

# output ------------------------------------------------------------------

write.csv(pkdat , file.path(outDataDir, "pkdat.csv") , row.names=FALSE, quote=FALSE)
write.csv(covdat, file.path(outDataDir, "covdat.csv"), row.names=FALSE, quote=FALSE)


