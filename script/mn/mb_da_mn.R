
# Library
library(tidyverse)
library(here)
library(readxl)
library(vegan)
library(stringr)

# Directory
sourceDataDir <- here("data/source/mn-biopsy-study")
dataDir <- here("data/derived/mn")

# Load metadata -----------------------------------------------------------

metadata <- read.csv(here(sourceDataDir, "CTSI_MN_master.csv"), header = TRUE) %>% 
  rename_all(tolower) %>%
  rename(subject_id = "patient.id") %>%
  # Missing WT: Subject 002- 87.5 kg, Subject 003- 91.4 kg
  mutate(
    weight = case_when(
      subject_id == 2 ~ 87.5,
      subject_id == 3 ~ 91.4,
      TRUE ~ weight)
  )

# Load microbiome genus data ----------------------------------------------

mb_genus <- read_xlsx(file.path(sourceDataDir, 
                                "HIV 16S rRNA Data (1).xlsx"), 
                      sheet="Genus", col_names=TRUE, range="A267:V531") %>% 
  # Rename UMN-021 to UMN-022 in the 16s dataset
  rename(MN_022="MN_021") %>% 
  rename(name="Percent") %>% 
  mutate(name=paste0("genus_", name)) %>% 
  # Keep xxxxxaceae_unclassified, remove other _unclassified
  filter(!grepl("_unclassified$", name) | grepl("aceae_unclassified$", name)) %>%
  #mutate(name = ifelse(grepl("_unclassified$", name), 
                       #paste0(sub("_unclassified$", "", name), "(f)"),
                       #name)) %>%
  column_to_rownames(var="name")

colnames(mb_genus) <- gsub("^MN_", "", colnames(mb_genus))

mb_genus2 <- mb_genus %>% data.table::transpose() # Transpose genus data

colnames(mb_genus2) <- rownames(mb_genus)
rownames(mb_genus2) <- colnames(mb_genus)

mb_genus <- mb_genus2 %>% 
  rownames_to_column(var = "subject_id")%>% 
  mutate(subject_id = as.integer(subject_id))
rm(mb_genus2)

names(mb_genus) <- tolower(names(mb_genus))


# Calculate microbiome alpha index ----------------------------------------------

otu_counts <- read_xlsx(file.path(sourceDataDir,
                                  "HIV 16S rRNA Data (1).xlsx"),
                        sheet="Taxonomy (norm=30000)", col_names=TRUE) %>% 
  # Rename UMN-021 to UMN-022 in the 16s dataset
  rename(MN_022="MN_021") %>% 
  select(starts_with("MN_")) 

colnames(otu_counts) <- str_remove(colnames(otu_counts), "^MN_")

shannon <- diversity(t(otu_counts), index = "shannon")

chao <- estimateR(t(otu_counts))["S.chao1", ]

mb_alpha <- data.frame(
  subject_id = as.numeric(colnames(otu_counts)),
  alpha_shannon = shannon,
  alpha_chao = as.numeric(chao)
)

# Load ftc drug exposures -------------------------------------------------

ftc_expo <- read.csv(file.path(dataDir, "ftc/ftc_expo.csv"))
ftc_etas <- read.csv(file.path(dataDir, "ftc/ftc_etas.csv"))

ftc <- ftc_etas %>% 
  dplyr::select(ID, subject_id=SUBJECT_ID) %>% # Get `subject_id`
  left_join(ftc_expo) %>% 
  rename_all(tolower) %>% 
  rename_at(vars(ka:frac24_2_d2), ~paste0("ftc_",.x))

# Merge data --------------------------------------------------------------

df <- left_join(metadata, ftc) %>% relocate(id, .after = subject_id)
df <- left_join(df, mb_alpha)
df <- left_join(df, mb_genus)
names(df)

# Output ------------------------------------------------------------------

write.csv(df, here(dataDir, "data.csv"), row.names = FALSE)
