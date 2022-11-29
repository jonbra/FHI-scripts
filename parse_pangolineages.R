library(tidyverse)

# Create list of BA.5 abbreviations for parsing in R ----------------------


# Read lineage descriptions from GitHub
pango <- read_delim(file = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")

# Extract only BA.5
pango %>% 
  # Get the BA.5's
  filter(str_detect(Description, "B.1.1.529.5")) %>% 
  # Drop BA.5 for now
  filter(Lineage != "BA.5") %>% 
  # Get the abbreviations
  mutate(tmp = str_sub(Lineage, 1, 3)) %>% 
  distinct(tmp) %>% 
  arrange(tmp) %>% print(n = 34)


# Create list of BA.5 lineages for the Nextstrain build file
pango_str <- pango %>% 
  # Get the BA.5's
  filter(str_detect(Description, "B.1.1.529.5")) %>% 
  filter(str_detect(Lineage, "\\*", negate = TRUE)) %>% 
  pull(Lineage)

# This almost does it. Copy and paste and fix the ends manually  
str_c("', '", pango_str, collapse = "")

