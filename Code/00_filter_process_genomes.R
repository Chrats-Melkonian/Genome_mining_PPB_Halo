# Set working directories
setwd(file.path(Sys.getenv("HOME"), "Desktop/Git/Genome_mining_PPB_Halo"))
data_all_spp <- read.csv("Data/data_all_spp.csv")

# Check unique species names
unique(unlist(lapply(strsplit(data_all_spp$Species, " "), function(x) x[1])))

# Species notes
# Rhabdochromatium: only 1 genome
# Isochromatium: no genome; no NCBI
# Ectothiorhodospira: no genome; found in NCBI

# Load necessary libraries
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(janitor)

# Set working directory for metadata files
setwd(paste(getwd(),"/Data/metadata/",sep = ""))
temp <- list.files()
table <- c()

# Read metadata files and combine them into a table
for (i in 1:length(temp)) {
  temp1 <- read.csv(temp[i])
  table <- rbind(table, temp1)
}

# Function to filter complete genomes based on quality and number of contigs
filtered_data <- function(data_name) {
  data <- data_name %>%
    clean_names() %>%
    mutate(genome_id = as.character(genome_id)) %>%
    filter((genome_status == 'Complete') %>% replace_na(TRUE)) %>%
    filter((genome_quality == 'Good') %>% replace_na(TRUE)) %>%
    mutate(plasmids = replace_na(plasmids, 0)) %>%
    mutate(nm_contigs = contigs - plasmids) %>%
    filter(nm_contigs <= 100)
  
  nm_sp <- nrow(data)
  return(data.frame(data))
}

filtered_complete1 <- filtered_data(table)

# Function to filter WGS genomes based on completeness and contamination
wgs_data <- function(data_name, i) {
  data <- data_name %>%
    clean_names() %>%
    filter((genome_status == "WGS") %>% replace_na(TRUE)) %>%
    filter((check_m_completeness > 90) %>% replace_na(TRUE)) %>%
    filter((check_m_contamination < 10) %>% replace_na(TRUE)) %>%
    mutate(plasmids = replace_na(plasmids, 0)) %>%
    mutate(nm_contigs = contigs - plasmids) %>%
    filter(nm_contigs <= 300)
  
  nm_sp <- nrow(data)
  return(data.frame(data))
}

# Filter WGS genomes
filtered_wgs <- wgs_data(table)

# Merge filtered complete genomes and WGS genomes
data_w <- merge(filtered_complete1, filtered_wgs, by = intersect(names(filtered_complete1), names(filtered_wgs)), all = TRUE)

# Write results to CSV and text file
write.csv(data_w, file = "metadata_final_selection.csv")

# Script for downloading genomes using wget
# for i in `cat genome_list_alln.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.faa"; done

# Set working directory for genomic data
setwd("PATHtoGENOMES")
output_folder <- "PATHtoOUTPUT"
file_names <- list.files()

# Generate Prodigal commands for each genome file
cmd <- c()

write_command <- function(i) {
  ioa <- sub(".fna", "_aa.faa", i)
  ion <- sub(".fna", "_dna.fna", i)
  x <- paste("prodigal -i ", paste(getwd(), i, sep = "/"), " -a ", paste(output_folder, "aa", ioa, sep = "/"), 
             " -d ", paste(output_folder, "dna", ion, sep = "/"), " -m -q", sep = "")
  return(x)
}

# Run Prodigal for each genome file
for (i in file_names) {
  cmd <- c(cmd, write_command(i))
}

sapply(cmd, system)

# Function to create EGGNOG command for annotation
cmdCreate <- function(input, output) {
  paste("emapper.py -i ", input, " -o ", output, " --cpu 45 --itype proteins --data_dir /home/chrats/NOBINFBACKUP/EGGNOG --pident 50 --query_cover 50 --subject_cover 50", sep = "")
}

# Set input and output directories for EGGNOG
output <- "PATHtoOUTPUTeggnog"
input <- "PATHtoAA"
setwd(input)
cmds <- c()
flist <- list.files(pattern = ".faa")

# Generate EGGNOG commands for each file
for (i in flist) {
  j <- paste(input, "/", i, sep = "")
  o <- paste(output, "/", gsub("_aa.faa|.faa", "_eggnog", i), sep = "")
  cmds <- c(cmds, cmdCreate(j, o))
}

# Run EGGNOG commands
sapply(cmds, system)
