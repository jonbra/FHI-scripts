pacman::p_load(tidyverse, Rsamtools, lubridate)

# Open the fasta index connection from the script index_fasta.R
# There needs to be a file with the same filename but ending with .fai in the folder
fa <- FaFile("/media/jonr/SATA6TB1/Gisaid/sequences.fasta")

# Make a GRanges object
gr <- as(seqinfo(fa), "GRanges")

# Load the metadata from the Entire Gisaid database
metadata_Gisaid <- read_tsv("/media/jonr/SATA6TB1/Gisaid/metadata.tsv")

# Sett minimumsdato
# Denne trekker 2 måneder fra dagens dato
#min_date <- Sys.Date() %m-% months(2)

# Filter the metadata
metadata_filtered <- metadata_Gisaid %>% 
  # Keep omicron only
  # filter(str_detect(`Pango lineage`, "BA.*") | str_detect(`Pango lineage`, "B.1.1.529")) %>% 
  # Drop BA.1
  # filter(str_detect(`Pango lineage`, "BA.1.*", negate = TRUE))
  # Keep BA.5 and BA.2.75*
  filter(str_detect(`Pango lineage`, "^BA.2.75.2.*") | str_detect(`Pango lineage`, "^BA.5.*") | str_detect(`Pango lineage`, "^BE.*") | str_detect(`Pango lineage`, "^BK.*") | str_detect(`Pango lineage`, "^BF.*") | str_detect(`Pango lineage`, "^BQ.*"))

# Clean up
rm(metadata_Gisaid)

# Match the names of the metadata with the fasta headers
new_gr <- getSeq(fa, gr[which(gsub("\\|.*", "", names(gr)) %in% metadata_filtered$`Virus name`)])

# Write files
outfile <- paste0("/media/jonr/SATA6TB1/Gisaid/", Sys.Date(), "_Gisaid.fasta")
# Write the fasta file
writeXStringSet(new_gr, outfile, format = "fasta")
# Write metadata
write_tsv(metadata_filtered, paste0("/media/jonr/SATA6TB1/Gisaid/", Sys.Date(), "_Gisaid.metadata.tsv"))

# Now, go to the script "get_data_from_BN.R" without removing any objects or variables from the environment
