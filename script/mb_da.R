
# Library
library(tidyverse)
library(here)
library(readxl)

# Directory
sourceDataDir <- here("data/source")
dataDir <- here("data/derived")

# Load metadata -----------------------------------------------------------

metadata <- read_xlsx(file.path(sourceDataDir, 
                                "Vaginal Microbiome Data Analysis.xlsx"), 
                      sheet="Metadata", col_names=TRUE) %>% 
  rename_all(tolower)

# Load microbiome phylum data ---------------------------------------------

mb_phylum <- read_xlsx(file.path(sourceDataDir, 
                                 "Vaginal Microbiome Data Analysis.xlsx"), 
                       sheet="Phylum", col_names=TRUE, range="A22:AY40") %>% 
  rename(name=`...1`) %>% 
  mutate(name=paste0("phylum_", name)) %>% 
  column_to_rownames(var="name") 

mb_phylum2 <- mb_phylum %>% data.table::transpose() # Transpose phylum data

colnames(mb_phylum2) <- rownames(mb_phylum)
rownames(mb_phylum2) <- colnames(mb_phylum)

mb_phylum <- mb_phylum2 %>% rownames_to_column(var = "msi_id")
rm(mb_phylum2)

names(mb_phylum) <- tolower(names(mb_phylum))

# Load microbiome genus data ----------------------------------------------

mb_genus <- read_xlsx(file.path(sourceDataDir, 
                                "Vaginal Microbiome Data Analysis.xlsx"), 
                      sheet="Genus", col_names=TRUE, range="A185:AY366") %>% 
  rename(name=`...1`) %>% 
  mutate(name=paste0("genus_", name)) %>% 
  column_to_rownames(var="name") 

mb_genus2 <- mb_genus %>% data.table::transpose() # Transpose genus data

colnames(mb_genus2) <- rownames(mb_genus)
rownames(mb_genus2) <- colnames(mb_genus)

mb_genus <- mb_genus2 %>% rownames_to_column(var = "msi_id")
rm(mb_genus2)

names(mb_genus) <- tolower(names(mb_genus))

# Load microbiome alpha data ----------------------------------------------

mb_alpha <- read_xlsx(file.path(sourceDataDir, 
                             "Vaginal Microbiome Data Analysis.xlsx"), 
                   sheet="Alpha", col_names=TRUE, range="A1:G51") %>% 
  rename(msi_id=group) %>% 
  dplyr::select(msi_id, shannon, chao) %>% 
  rename_at(vars(shannon:chao), ~paste0("alpha_", .x)) %>% 
  rename_all(tolower)

# Load tfv drug exposures -------------------------------------------------

tfv_expo <- read.csv(file.path(dataDir, "tfv/tfv_expo.csv"))
tfv_etas <- read.csv(file.path(dataDir, "tfv/tfv_etas.csv"))

tfv <- tfv_etas %>% 
  dplyr::select(ID, subject_id=SUBJECT_ID) %>% # Get `subject_id`
  left_join(tfv_expo) %>% 
  rename_all(tolower) %>% 
  rename_at(vars(ka:frac24_3_d2), ~paste0("tfv_",.x))

# Load 3tc drug exposures -------------------------------------------------

ttc_expo <- read.csv(file.path(dataDir, "3tc/3tc_expo.csv"))
ttc_etas <- read.csv(file.path(dataDir, "3tc/3tc_etas.csv"))

ttc <- ttc_etas %>% 
  dplyr::select(ID, subject_id=SUBJECT_ID) %>% # Get `subject_id`
  left_join(ttc_expo) %>% 
  rename_all(tolower) %>% 
  rename_at(vars(ka:cave_t_d2), ~paste0("3tc_",.x))

# Merge data --------------------------------------------------------------

df <- left_join(metadata, tfv) %>% left_join(ttc) %>% dplyr::select(-id)
df <- left_join(df, mb_alpha)
df <- left_join(df, mb_phylum)
df <- left_join(df, mb_genus)
names(df)

# Output ------------------------------------------------------------------

write.csv(df, here(dataDir, "data.csv"), row.names = FALSE)
