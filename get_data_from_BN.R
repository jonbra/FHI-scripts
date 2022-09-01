pacman::p_load(tidyverse, phylotools)

# Create empty objects for the fasta sequences
FHI_fastas <- tibble(
  "seq.name" = character(),
  "seq.text" = character()
)
MIK_fastas <- tibble(
  "seq.name" = character(),
  "seq.text" = character()
)
Nano_fastas <- tibble(
  "seq.name" = character(),
  "seq.text" = character()
)
Artic_fastas <- tibble(
  "seq.name" = character(),
  "seq.text" = character()
)

# Load BN (remember to refresh first)
load("/mnt/N/Virologi/JonBrate/Prosjekter/BN.RData")

# Convert empty strings to NA and clean up
BN <- BN %>% mutate_all(list(~na_if(.,""))) %>% 
  # Endre Trøndelag til Trondelag
  mutate("FYLKENAVN" = str_replace(FYLKENAVN, "Tr\xf8ndelag", "Trondelag")) %>%
  # Endre Møre og Romsdal
  mutate("FYLKENAVN" = str_replace(FYLKENAVN, "M\xf8re", "More")) %>%
  # Endre Sør
  mutate("FYLKENAVN" = str_replace(FYLKENAVN, "S\xf8r", "Sor")) %>% 
  # Fix date format
  mutate("PROVE_TATT" = lubridate::ymd(PROVE_TATT)) %>%
  # Drop samples witout collection date
  filter(!is.na(PROVE_TATT))

# Sett Virus name som fasta header
# Først lage en mapping mellom KEY og virus name
SEQUENCEID_virus_mapping_FHI <- BN %>%
  filter(PROVE_TATT >= "2022-01-01") %>% 
  # Keep BA.5 only
  filter(str_detect(PANGOLIN_NOM, "^BA.5.*")) %>% 
  # Keep omicron only
  #filter(str_detect(PANGOLIN_NOM, "BA.*") | str_detect(PANGOLIN_NOM, "B.1.1.529")) %>% 
  # Drop BA.1
  #filter(str_detect(PANGOLIN_NOM, "BA.1.*", negate = TRUE)) %>% 
  filter(is.na(GISAID_EPI_ISL)) %>% 
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "FHI")) %>% 
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  # Rename mis-spelled plates
  mutate(SEKV_OPPSETT_SWIFT7 = str_replace(SEKV_OPPSETT_SWIFT7, "FHI432s", "FHI432-S")) %>% 
  mutate(SEKV_OPPSETT_SWIFT7 = str_replace(SEKV_OPPSETT_SWIFT7, "FHI429n", "FHI429-N")) %>% 
  select(KEY, SEQUENCEID_SWIFT, covv_virus_name, SEKV_OPPSETT_SWIFT7)

SEQUENCEID_virus_mapping_MIK <- BN %>%
  filter(PROVE_TATT >= "2022-01-01") %>% 
  # Keep BA.5 only
  filter(str_detect(PANGOLIN_NOM, "^BA.5.*")) %>% 
  # Keep omicron only
  #filter(str_detect(PANGOLIN_NOM, "BA.*") | str_detect(PANGOLIN_NOM, "B.1.1.529")) %>% 
  # Drop BA.1
  #filter(str_detect(PANGOLIN_NOM, "BA.1.*", negate = TRUE)) %>% 
  filter(is.na(GISAID_EPI_ISL)) %>% 
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>% 
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1)) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(SEQUENCEID_SWIFT, KEY, covv_virus_name, SEKV_OPPSETT_SWIFT7)

SEQUENCEID_virus_mapping_Artic <- BN %>%
  filter(PROVE_TATT >= "2022-01-01") %>% 
  # Keep BA.5 only
  filter(str_detect(PANGOLIN_NOM, "^BA.5.*")) %>% 
  # Keep omicron only
  #filter(str_detect(PANGOLIN_NOM, "BA.*") | str_detect(PANGOLIN_NOM, "B.1.1.529")) %>% 
  # Drop BA.1
  #filter(str_detect(PANGOLIN_NOM, "BA.1.*", negate = TRUE)) %>% 
  filter(is.na(GISAID_EPI_ISL)) %>% 
  filter(str_detect(RES_CDC_INFB_CT, "Artic")) %>%
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(KEY, RES_CDC_INFB_CT, covv_virus_name, SAMPLE_CATEGORY)

SEQUENCEID_virus_mapping_Nano <- BN %>%
  filter(PROVE_TATT >= "2022-01-01") %>% 
  # Keep BA.5 only
  filter(str_detect(PANGOLIN_NOM, "^BA.5.*")) %>% 
  # Keep omicron only
  #filter(str_detect(PANGOLIN_NOM, "BA.*") | str_detect(PANGOLIN_NOM, "B.1.1.529")) %>% 
  # Drop BA.1
  #filter(str_detect(PANGOLIN_NOM, "BA.1.*", negate = TRUE)) %>% 
  filter(is.na(GISAID_EPI_ISL)) %>% 
  filter(str_detect(SEKV_OPPSETT_NANOPORE, "Nano")) %>%
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(KEY, SEQUENCEID_NANO29, covv_virus_name, SEKV_OPPSETT_NANOPORE)

# Get the different "oppsett" to search only specific folders for sequences
FHI_plater <- SEQUENCEID_virus_mapping_FHI %>% 
  select(SEKV_OPPSETT_SWIFT7, SEQUENCEID_SWIFT) %>% 
  drop_na(SEKV_OPPSETT_SWIFT7) %>% 
  distinct(SEKV_OPPSETT_SWIFT7) %>% 
  pull(SEKV_OPPSETT_SWIFT7)

FHI_IDs <- SEQUENCEID_virus_mapping_FHI %>% 
  select(SEKV_OPPSETT_SWIFT7, SEQUENCEID_SWIFT) %>% 
  drop_na(SEKV_OPPSETT_SWIFT7) %>% 
  filter(str_detect(SEQUENCEID_SWIFT, "konklu", negate = TRUE))

MIK_plater <- SEQUENCEID_virus_mapping_MIK %>% 
  select(SEKV_OPPSETT_SWIFT7, SEQUENCEID_SWIFT) %>% 
  drop_na(SEKV_OPPSETT_SWIFT7) %>% 
  distinct(SEKV_OPPSETT_SWIFT7) %>% 
  pull(SEKV_OPPSETT_SWIFT7)

MIK_IDs <- SEQUENCEID_virus_mapping_MIK %>% 
  select(SEKV_OPPSETT_SWIFT7, SEQUENCEID_SWIFT) %>% 
  drop_na(SEKV_OPPSETT_SWIFT7) %>% 
  filter(str_detect(SEQUENCEID_SWIFT, "konklu", negate = TRUE))

Artic_plater <- SEQUENCEID_virus_mapping_Artic %>% 
  select(RES_CDC_INFB_CT, SAMPLE_CATEGORY) %>% 
  drop_na(SAMPLE_CATEGORY) %>% 
  separate(SAMPLE_CATEGORY, into = c(NA, "oppsett_nr", NA), sep = "/") %>% 
  distinct(oppsett_nr) %>% 
  pull(oppsett_nr)

Artic_IDs <- SEQUENCEID_virus_mapping_Artic %>% 
  select(RES_CDC_INFB_CT, SAMPLE_CATEGORY) %>% 
  drop_na(SAMPLE_CATEGORY) %>% 
  filter(str_detect(SAMPLE_CATEGORY, "konklu", negate = TRUE))

Nano_plater <- SEQUENCEID_virus_mapping_Nano %>% 
  select(SEKV_OPPSETT_NANOPORE, SEQUENCEID_NANO29) %>% 
  drop_na(SEKV_OPPSETT_NANOPORE) %>% 
  separate(SEKV_OPPSETT_NANOPORE, into = c("oppsett_nr", NA), sep = "/") %>% 
  mutate(oppsett_nr = str_remove(oppsett_nr, "Nr")) %>% 
  distinct(oppsett_nr) %>% 
  pull(oppsett_nr)

Nano_IDs <- SEQUENCEID_virus_mapping_Nano %>% 
  select(SEKV_OPPSETT_NANOPORE, SEQUENCEID_NANO29) %>% 
  drop_na(SEQUENCEID_NANO29) %>% 
  filter(str_detect(SEQUENCEID_NANO29, "konklu", negate = TRUE))


#list directories to search for sequences
dirs_fhi <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2021/", recursive = FALSE), list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2022/", recursive = FALSE))
dirs_mik <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK", recursive = FALSE)
dirs_artic <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2021", recursive = FALSE), list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2022", recursive = FALSE))
dirs_nano <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2021", recursive = FALSE), list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2022", recursive = FALSE))

# Get the sequences
# SWIFT FHI
for (i in seq_along(FHI_plater)) {
  # Pick our the relevant oppsett
  dir <- dirs_fhi[grep(FHI_plater[i], dirs_fhi)]
  
  # List the files
  filepaths <- list.files(path = dir,
                          pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$",
                          full.names = TRUE,
                          recursive = TRUE)
  # Clean ut filepaths for grepping later (in cases when "_2")
  names <- str_remove(filepaths, "_ivar\\.consensus\\.masked_Nremoved\\.fa$")
  
  # Find which filepaths to keep
  seq_ids <- FHI_IDs %>% 
    filter(str_detect(SEKV_OPPSETT_SWIFT7, FHI_plater[i]))
  
  keep <- vector("character")
  for (y in seq_along(seq_ids$SEQUENCEID_SWIFT)){
    try(keep[y] <- filepaths[grep(paste0(seq_ids$SEQUENCEID_SWIFT[y], "\\b"), names)]) 
  }
  # Drop empty elements in keep (where no sequence file was found for the sequence id)
  keep <- keep[!is.na(keep)]
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  if (length(keep) > 0) {
    for (f in seq_along(keep)){
      tmp <- read.fasta(keep[f])      # read the file
      fastas <- rbind(fastas, tmp)    # append the current file
    }
    # Convert to tibble for easier manipulation
    fastas <- as_tibble(fastas)
    
    # Fix names to match SEQUENCEID_SWIFT
    fastas <- fastas %>%
      mutate(SEQUENCEID_SWIFT = str_remove(seq.name, "_ivar_masked"))
    
    # Add virus name
    fastas <- left_join(fastas, SEQUENCEID_virus_mapping_FHI, by = "SEQUENCEID_SWIFT") %>%
      select(`seq.name` = covv_virus_name,
             seq.text)
    # Join with final object
    FHI_fastas <- bind_rows(FHI_fastas, fastas)
  }
}
# SWIFT MIK
for (i in seq_along(MIK_plater)) {
  # Pick our the relevant oppsett
  dir <- dirs_mik[grep(MIK_plater[i], dirs_mik)]
  
  # List the files
  filepaths <- list.files(path = dir,
                          pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$",
                          full.names = TRUE,
                          recursive = TRUE)
  # Clean ut filepaths for grepping later (in cases when "_2")
  names <- str_remove(filepaths, "_ivar\\.consensus\\.masked_Nremoved\\.fa$")
  
  # Find which filepaths to keep
  seq_ids <- MIK_IDs %>% 
    filter(str_detect(SEKV_OPPSETT_SWIFT7, MIK_plater[i]))
  
  keep <- vector("character")
  for (y in seq_along(seq_ids$SEQUENCEID_SWIFT)){
    try(keep[y] <- filepaths[grep(paste0(seq_ids$SEQUENCEID_SWIFT[y], "\\b"), names)]) 
  }
  # Drop empty elements in keep (where no sequence file was found for the sequence id)
  keep <- keep[!is.na(keep)]
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match SEQUENCEID_SWIFT
  fastas <- fastas %>%
    mutate(SEQUENCEID_SWIFT = str_remove(seq.name, "_ivar_masked"))
  
  # Add virus name
  fastas <- left_join(fastas, SEQUENCEID_virus_mapping_MIK, by = "SEQUENCEID_SWIFT") %>%
    select(`seq.name` = covv_virus_name,
           seq.text)
  # Join with final object
  MIK_fastas <- bind_rows(MIK_fastas, fastas)
}
# Artic Illumina
for (i in seq_along(Artic_plater)) {
  # Pick our the relevant oppsett
  dir <- dirs_artic[grep(Artic_plater[i], dirs_artic)]
  
  # List the files
  filepaths <- list.files(path = dir,
                          pattern = "consensus\\.fa$",
                          full.names = TRUE,
                          recursive = TRUE)
  
  # Find which filepaths to keep
  seq_ids <- Artic_IDs %>% 
    filter(str_detect(SAMPLE_CATEGORY, Artic_plater[i]))
  
  keep <- vector("character")
  for (y in seq_along(seq_ids$RES_CDC_INFB_CT)){
    try(keep[y] <- filepaths[grep(seq_ids$RES_CDC_INFB_CT[y], filepaths)])
  }
  # Drop empty elements in keep (where no sequence file was found for the sequence id)
  keep <- keep[!is.na(keep)]
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  if (length(keep) > 0) {
    for (f in seq_along(keep)){
      tmp <- read.fasta(keep[f])      # read the file
      fastas <- rbind(fastas, tmp)    # append the current file
    }
    # Convert to tibble for easier manipulation
    fastas <- as_tibble(fastas)
    
    # Add virus name
    fastas <- left_join(fastas, SEQUENCEID_virus_mapping_Artic, by = c("seq.name" = "RES_CDC_INFB_CT")) %>%
      select(`seq.name` = covv_virus_name,
             seq.text)
    # Join with final object
    Artic_fastas <- bind_rows(Artic_fastas, fastas)
  }
}
# Artic Nanopore
for (i in seq_along(Nano_plater)) {
  # Pick our the relevant oppsett
  dir <- dirs_nano[grep(Nano_plater[i], dirs_nano)]
  
  # List the files
  filepaths <- list.files(path = dir,
                          pattern = "consensus\\.fasta$",
                          full.names = TRUE,
                          recursive = TRUE)
  
  # Find which filepaths to keep
  seq_ids <- Nano_IDs %>% 
    filter(str_detect(SEKV_OPPSETT_NANOPORE, Nano_plater[i]))
  
  keep <- vector("character")
  for (y in seq_along(seq_ids$SEQUENCEID_NANO29)){
    try(keep[y] <- filepaths[grep(seq_ids$SEQUENCEID_NANO29[y], filepaths)])
  }
  # Drop empty elements in keep (where no sequence file was found for the sequence id)
  keep <- keep[!is.na(keep)]
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match SEQUENCEID_NANO29
  fastas <- fastas %>%
    separate("seq.name", into = c("SEQUENCEID_NANO29", NA, NA), sep = "/", remove = F)
  
  # Add virus name
  fastas <- left_join(fastas, SEQUENCEID_virus_mapping_Nano, by = "SEQUENCEID_NANO29") %>%
    select(`seq.name` = covv_virus_name,
           seq.text)
  # Join with final object
  Nano_fastas <- bind_rows(Nano_fastas, fastas)
}

# Remove excact duplicates 
FHI_fastas <- distinct(FHI_fastas)
MIK_fastas <- distinct(MIK_fastas)
Artic_fastas <- distinct(Artic_fastas)
Nano_fastas <- distinct(Nano_fastas)

# Join all
Total_fastas <- bind_rows(
  FHI_fastas,
  MIK_fastas,
  Artic_fastas,
  Nano_fastas
)

SEQUENCEID_virus_mapping <- bind_rows(
  SEQUENCEID_virus_mapping_FHI,
  SEQUENCEID_virus_mapping_MIK,
  SEQUENCEID_virus_mapping_Artic,
  SEQUENCEID_virus_mapping_Nano
)


# First get the Eksterne metadata
eksterne_meta <- BN %>%
  filter(PROVE_TATT >= "2022-01-01") %>% 
  # Keep BA.5 only
  filter(str_detect(PANGOLIN_NOM, "^BA.5.*")) %>% 
  # Keep omicron only
  #filter(str_detect(PANGOLIN_NOM, "BA.*") | str_detect(PANGOLIN_NOM, "B.1.1.529")) %>% 
  # Drop BA.1
  #filter(str_detect(PANGOLIN_NOM, "BA.1.*", negate = TRUE)) %>%  
  filter(is.na(GISAID_EPI_ISL)) %>% 
  filter(str_detect(KEY, "SUS") | str_detect(KEY, "STO") | str_detect(KEY, "UNN") | str_detect(KEY, "HUS")) %>% 
  # Select final columns
  select("strain" = KEY,
         "date" = PROVE_TATT,
         "division" = FYLKENAVN,
         PANGOLIN_NOM) %>% 
  add_column("region" = "Europe",
             "country" = "Norway")

# Then get SUS fastas
SUS_files <- list.files("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/SUS/",
                        pattern = ".fasta$",
                        full.names = TRUE,
                        recursive = TRUE)

# Get the timestamps
file_info <- file.info(SUS_files, extra_cols = FALSE)
files_to_keep <- as_tibble(rownames_to_column(file_info, var = "filename")) %>% 
  filter(mtime >= "2022-01-01") %>% 
  pull(filename)

# First create empty data frame to fill
SUS_fastas <- data.frame(seq.name = character(),
                         seq.text = character())

for (f in seq_along(files_to_keep)){
  tmp <- read.fasta(files_to_keep[f])      # read the file
  SUS_fastas <- rbind(SUS_fastas, tmp)    # append the current file
}

SUS_fastas <- as_tibble(SUS_fastas)

SUS_fastas <- SUS_fastas %>% 
  mutate(tmp = str_extract(seq.name, "H2.*")) %>% 
  filter(!is.na(tmp)) %>% 
  mutate(tmp2 = str_remove(tmp, "\\.consensus\\.fasta")) %>% 
  mutate(tmp3 = paste0("SUS-", tmp2)) 

SUS_fastas <- left_join(eksterne_meta, SUS_fastas, by = c("strain" = "tmp3")) %>% 
  select("seq.name" = strain,
         "seq.text") 
SUS_fastas <- SUS_fastas %>% drop_na(seq.text)

# Get the HUS fastas
HUS_files <- list.files("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/HUS/",
                        pattern = ".fasta$",
                        full.names = TRUE,
                        recursive = TRUE)

# First create empty data frame to fill
HUS_fastas <- data.frame(seq.name = character(),
                         seq.text = character())

for (f in seq_along(HUS_files)){
  tmp <- read.fasta(HUS_files[f])      # read the file
  HUS_fastas <- rbind(HUS_fastas, tmp)    # append the current file
}

HUS_fastas <- as_tibble(HUS_fastas)

HUS_fastas <- left_join(eksterne_meta, HUS_fastas, by = c("strain" = "seq.name")) %>% 
  select("seq.name" = strain,
         "seq.text") %>% 
  drop_na(seq.text)

# Merge with Total fastas
Total_fastas <- bind_rows(Total_fastas,
                          SUS_fastas,
                          HUS_fastas)


# Get metadata from BN
metadata_BN <- left_join(SEQUENCEID_virus_mapping, BN, by = "KEY") %>% 
  # Select final columns
  select("strain" = covv_virus_name,
          "date" = PROVE_TATT,
         "division" = FYLKENAVN,
         PANGOLIN_NOM) %>% 
  add_column("region" = "Europe",
             "country" = "Norway")

# Join with eksterne
metadata_BN <- bind_rows(
  metadata_BN,
  eksterne_meta
) %>% 
  # Adding dummy date submitted
  add_column("date_submitted" = Sys.Date())

# Create the same columns as the Gisaid metadata used - just in case
# NB! Må legge inn disse kolonene? https://github.com/nextstrain/augur/issues/575
# country division location region_exposure country_exposure division_exposure segment
# og country_confidence?
metadata_BN <- metadata_BN %>%
  unite("Location", c(region, country, division), sep = " / ", remove = FALSE) %>% 
  add_column("Type" = "betacoronavirus",
             "Accession ID" = NA,
             "Additional location information" = "unknown",
             "Sequence length" = NA,
             "Host" = "Human",
             "Patient age" = "unknown",
             "Gender" = "unknown",
             "Clade" = "unknown",
             "Pangolin version" = "unknown",
             "Variant" = "unknown",
             "AA Substitutions" = "unknown",
             "Is reference?" = NA,
             "Is complete?" = NA,
             "Is high coverage?" = "True",
             "Is low coverage?" = "False",
             "N-Content" = "unknown",
             "GC-Content" = "unknown") %>% 
  select("Virus name" = strain,
         Type,
         "Accession ID",
         "Collection date" = date,
         Location,
         "Additional location information",
         "Sequence length",
         Host,
         "Patient age",
         Gender,
         Clade,
         "Pango lineage" = PANGOLIN_NOM,
         "Pangolin version",
         "Variant",
         "AA Substitutions",
         "Submission date" = date_submitted,
         "Is reference?",
         "Is complete?",
         "Is high coverage?",
         "Is low coverage?",
         "N-Content",
         "GC-Content")
metadata_BN <- metadata_BN %>% 
  mutate(`Collection date` = as.character(`Collection date`)) 


#Gisaid_metadata <- metadata_filtered %>% 
#  # Massage columns
#  mutate(`Pangolin version` = as.character(`Pangolin version`),
#         `Sequence length` = as.character(`Sequence length`),
#         `Submission date` = as.character(`Submission date`),
#         `Is high coverage?` = as.character(`Is high coverage?`),
#         `Is low coverage?` = as.character(`Is low coverage?`),
#         `N-Content` = as.character(`N-Content`),
#         `GC-Content` = as.character(`GC-Content`))

# Check if metadata BN is identical to metadata Gisaid
# Merge with metadata_filtered from parse_large_fasta.R

# Write files
#if (identical(colnames(Gisaid_metadata), colnames(metadata_BN))) {
  outfile_fasta <- paste0("/home/jonr/Prosjekter/Nextstrain_mamba/ncov/data/", 
                          (format(Sys.Date(), "%Y.%m.%d")),
                          "/",
                          Sys.Date(),
                          "_BN.fasta")
  outfile_metadata <- paste0("/home/jonr/Prosjekter/Nextstrain_mamba/ncov/data/", 
                             (format(Sys.Date(), "%Y.%m.%d")),
                             "/",
                             Sys.Date(),
                             "_BN.metadata.tsv")
  
  dat2fasta(Total_fastas, outfile = outfile_fasta)
  write_tsv(metadata_BN, file = outfile_metadata)
#}





# Make overlay file -------------------------------------------------------

# Get FHI data from BN from the last two months, including what's submitted to Gisaid
# Sett Virus name som fasta header
# Først lage en mapping mellom KEY og virus name
SEQUENCEID_virus_mapping_FHI <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "FHI")) %>% 
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(KEY, SEQUENCEID_SWIFT, covv_virus_name)

SEQUENCEID_virus_mapping_MIK <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>% 
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1)) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(SEQUENCEID_SWIFT, KEY, covv_virus_name)

SEQUENCEID_virus_mapping_Artic <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(RES_CDC_INFB_CT, "Artic")) %>%
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(KEY, RES_CDC_INFB_CT, covv_virus_name)

SEQUENCEID_virus_mapping_Nano <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(SEKV_OPPSETT_NANOPORE, "Nano")) %>%
  # Trenger også å lage Virus name
  # Lage kolonne for "year"
  separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
  # Trekke ut sifrene fra 5 og til det siste fra BN KEY
  mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
  # Fjerne ledende nuller fra stammenavnet
  mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>%
  # Legge til kolonner med fast informasjon for å lage "Virus name" senere
  add_column("Separator" = "/",
             "GISAID_prefix" = "hCoV-19/",
             "Country" = "Norway/",
             "Continent" = "Europe/") %>%
  # Make "Virus name" column
  unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
  select(KEY, SEQUENCEID_NANO29, covv_virus_name)

SEQUENCEID_virus_mapping_Ahus <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(KEY, "Ahus")) %>%
  # create a covv_virus_name identical to KEY
  mutate("covv_virus_name" = KEY) %>% 
  select(KEY, covv_virus_name)

SEQUENCEID_virus_mapping_SUS <- BN %>%
  # Extract all data from last two months
  filter(PROVE_TATT >= min_date) %>% 
  # Remove all Delta
  filter(str_detect(PANGOLIN_NOM, "AY.*", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "B.1.617", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "Failed", negate = TRUE)) %>% 
  filter(str_detect(PANGOLIN_NOM, "konklu", negate = TRUE)) %>% 
  # Remove NA, na, Na, etc.
  #filter(str_detect(PANGOLIN_NOM, "(?i)NA", negate = TRUE)) %>% 
  filter(str_detect(KEY, "SUS")) %>%
  # create a covv_virus_name identical to KEY
  mutate("covv_virus_name" = KEY) %>% 
  select(KEY, covv_virus_name)


# Join everything
SEQUENCEID_virus_mapping <- bind_rows(
  SEQUENCEID_virus_mapping_FHI,
  SEQUENCEID_virus_mapping_MIK,
  SEQUENCEID_virus_mapping_Artic,
  SEQUENCEID_virus_mapping_Nano,
  SEQUENCEID_virus_mapping_Ahus,
  SEQUENCEID_virus_mapping_SUS
)

# Strip hCoV-19
SEQUENCEID_virus_mapping <- SEQUENCEID_virus_mapping %>% 
  mutate("strain" = str_remove(covv_virus_name, "hCoV-19/")) 

# Then get metadata from BN
overlay <- left_join(SEQUENCEID_virus_mapping, BN, by = "KEY") %>% 
  # Select final columns
  select(strain,
         "date" = PROVE_TATT,
         "division" = FYLKENAVN,
         "Vaksine" = VAKSINE,
         P,
         "Pangolin_BN" = PANGOLIN_NOM,
         "Utbrudd" = SESONG,
         "Reise" = REISE,
         "Sted" = STED,
         "St" = ST,
         "Vaks" = VAKS,
         "Utenlands_Siste7dager" = UTENLANDS_SISTE7DAGER,
         "gisaid_epi_isl" = GISAID_EPI_ISL) %>% 
  add_column("region" = "Europe",
             "country" = "Norway")

overlay_filename <- paste0("/mnt/N/Virologi/JonBrate/Nextstrain/Omikron/auspice/", 
                        (format(Sys.Date(), "%Y.%m.%d")),
                        "-overlay.tsv")
write_tsv(overlay, file = overlay_filename)
