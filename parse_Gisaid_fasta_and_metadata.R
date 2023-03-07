pacman::p_load(tidyverse, Rsamtools, lubridate)

# Open the fasta index connection from the script index_fasta.R
# There needs to be a file with the same filename but ending with .fai in the folder
fa <- FaFile("/media/jonr/SATA6TB/Gisaid/sequences.fasta")

# Make a GRanges object
gr <- as(seqinfo(fa), "GRanges")

# Load the metadata from the Entire Gisaid database
metadata_Gisaid <- read_tsv("/media/jonr/SATA6TB/Gisaid/metadata.tsv")


# Read lineage descriptions from GitHub
pango <- read_delim(file = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")

# Create list of BA.5 and BA.2.75 lineages for the Nextstrain build file
pango_str <- pango %>% 
  # Get the BA.5's
  filter(str_detect(Description, "B.1.1.529.5") | str_detect(Description, "B.1.1.529.2.75")) %>% 
  # Remove some withdrawn lineages
  filter(str_detect(Lineage, "\\*", negate = TRUE)) %>% 
  # Pull all the aliases into a character vector
  pull(Lineage)

# Filter the metadata for BA.5 and BA.2
filtered <- metadata_Gisaid %>% 
  filter(`Pango lineage` %in% pango_str)

# Get the XBB.1.5 (this abbreviation is not linked to full pangos in the GitHub file. But used in the Gisaid metadata)
XBB <- metadata_Gisaid %>% 
  filter(str_detect(`Pango lineage`, "^XBB.1.5"))

# Combine the metadata
metadata_filtered <- bind_rows(filtered, XBB)

# Clean up
rm(metadata_Gisaid)
gc()

# Match the names of the metadata with the fasta headers
# Using AAStringSet to account for possible presence of Non-DNA letters. 
# The AAStringSet allows for all Letters
new_gr <- getSeq(fa, gr[which(gsub("\\|.*", "", names(gr)) %in% metadata_filtered$`Virus name`)], as="AAStringSet")

# Write files
outfile <- paste0("/home/jonr/Prosjekter/Nextstrain_mamba/ncov/data/SC2_weekly/", "Gisaid.fasta")
# Write the fasta file
writeXStringSet(new_gr, outfile, format = "fasta")
# Write metadata
write_tsv(metadata_filtered, paste0("/home/jonr/Prosjekter/Nextstrain_mamba/ncov/data/SC2_weekly/", "Gisaid.metadata.tsv"))

# Now, go to the script "get_data_from_BN.R" without removing any objects or variables from the environment
