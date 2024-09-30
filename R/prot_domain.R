
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
library(tools)
library(ggplot2)
library(tidyverse)
library(protr)
library(ragp)
library(yaml)

ini <- read_yaml("pipeline.yml")

folder_path <- ini$folder_path

### Prediction of protein domains 

output_dir <- "Protein_domains_isoform_ragp"
dir.create(output_dir, showWarnings = FALSE)

# Define the function to read sample names from each TSV file and perform domain prediction

prot_domain <- function(folder_path, file_name) {
  # Construct the full file path
  file_path <- file.path(folder_path, file_name)
  
  # Generate the plot file name
  plot_file <- file.path(output_dir, paste0(file_path_sans_ext(file_name), ".pdf"))
  
  # Skip processing if the output file already exists
  if (file.exists(plot_file)) {
    message(paste("Skipping", file_name, "as it has already been processed."))
    return(NULL)
  }
  
  # Read the TSV file
  data <- read.delim(file_path, sep = "\t", header = TRUE)
  
  # Filter data
  df_prot <- data %>%
    dplyr::filter(!is.na(uniprot_id))%>%
    dplyr::filter(!duplicated(uniprot_id))
  
  df_prot$uniprot_id <- gsub("\\..*", "",df_prot$uniprot_id)
  
  ids <- unique(df_prot$uniprot_id)
  
  seqs <- df_prot$protein_sequence 
  
  # Function to retry the plot_prot call
  retry_plot_prot <- function(seqs, ids, attempts = 5, timeout = 10) {
    p <- NULL
    for (i in 1:attempts) {
      try({
        p <- plot_prot(seqs, ids, hyp = FALSE, ag = FALSE, domain = "hmm")
        if (!is.null(p)) break
      }, silent = TRUE)
      Sys.sleep(timeout) # Wait before retrying
    }
    return(p)
  }
  
  # Call plot_prot with retries
  p <- retry_plot_prot(seqs, ids)
  
  # Check if plot creation was successful
  if (is.null(p)) {
    warning(paste("Failed to create plot for", file_name, "after multiple attempts. Skipping this file."))
    return(NULL)
  }
  
  ggsave(plot_file, plot = p, width = 10, height = 6)
}

# Apply the alignment function to each TSV file in the folder
# Example usage
# Define the folder path containing your TSV files
folder_path <- folder_path

# Get a list of TSV files in the folder
tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = FALSE)

lapply(tsv_files, function(file_name) prot_domain(folder_path, file_name))