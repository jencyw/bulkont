
### Author: Chen-Yi Wang
### Title: "Alignment and the domains of protein isoforms"
### Goal: to perform protein alignment and domain prediction of protein isoforms  
### Output: protein alignment and schematics of protein domains

### Before run, please check the following requirement:   
# DTU analysis has been done and the output files in the folder of StageR_df_plotproportions
# R availability in the environment
# Availability of the required R packages

### Run this pipeline in the folder of DTU_DRIMSeq.dir

library(dplyr)
library(stringr)
library(msa)
library(Biostrings)
library(tools)
library(ggplot2)
library(tidyverse)
library(yaml)

ini <- read_yaml("pipeline.yml")
folder_path <- ini$folder_path

### Prediction of protein domains

output_dir <- "Protein_AAAlignment"
dir.create(output_dir, showWarnings = FALSE)

# Define the function to read sample names from each TSV file and perform alignment
alignment <- function(folder_path, file_name) {
  # Construct the full file path
  file_path <- file.path(folder_path, file_name)
  
  # Generate the plot file name
  plot_file <- file.path(output_dir, paste0(file_path_sans_ext(file_name), "_AAAlignment", ".pdf"))
  
  # Skip processing if the output file already exists
  if (file.exists(plot_file)) {
    message(paste("Skipping", file_name, "as it has already been processed."))
    return(NULL)
  }
  
  # Read the TSV file
  data <- read.delim(file_path, sep = "\t", header = TRUE)
  
  # Filter data
  data <- data %>%
    filter(!is.na(protein_sequence)) %>%
    filter(!duplicated(tx_character))
  
  # Skip if filtered_data is empty
  if (nrow(data) <= 1) {
    warning(paste("data is either empty or only one type of protein sequence for file:", file_name, ". Skipping this file."))
    return(NULL)
  }
  
  # Extract and prepare sequences
  aa_seq <- data$protein_sequence
  names(aa_seq) <- data$tx_character
  
  # Create AAStringSet object
  aa_align <- AAStringSet(aa_seq, use.names = TRUE)
  
  # Perform multiple sequence alignment
  aa_multialign <- msa(aa_align, "ClustalOmega")
  
  # Save the alignment plot as a PDF with error handling
  tryCatch({
    msaPrettyPrint(aa_multialign, output = "pdf", showNames = "left", file = plot_file, showLogo = "none", askForOverwrite = FALSE, verbose = FALSE)
  }, error = function(e) {
    if (grepl("TeX capacity exceeded", e$message)) {
      warning(paste("TeX capacity exceeded for file:", file_name, ". Skipping this file."))
    } else {
      stop(e) # Re-throw the error if it's not a TeX capacity issue
    }
  })
}

# Apply the alignment function to each TSV file in the folder
# Get a list of TSV files in the folder
tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = FALSE)

# Apply the alignment function to each TSV file in the folder
lapply(tsv_files, function(file_name) alignment(folder_path, file_name))

# Define the source and destination directories
source_dir <- getwd()
dest_dir <- "Protein_AAAlignment"

# List all files in the source directory that contain "Alignment" in their names
files_to_move <- list.files(path = source_dir, pattern = "Alignment", full.names = TRUE)

# Move each file to the destination directory
for (file_path in files_to_move) {
  file_name <- basename(file_path)
  new_file_path <- file.path(dest_dir, file_name)
  
  # Move the file
  file.rename(from = file_path, to = new_file_path)
}

# Print the list of moved files
print(paste("Moved files:", files_to_move))
