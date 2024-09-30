
### Author: Chen-Yi Wang
### Title: "Transcript discovery using Bambu"
### Goal: to perform transcript identification and novel transcript discovery using multiple bam files  
### Output: count matrices including counts on gene and transcript level

### Before run, please check the following requirement:   
# R availability in the environment
# Availability of the required R packages: Bambu

library(devtools)
library(Rsamtools)
library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)
library(NanoporeRNASeq)
library(yaml)

ini <- read_yaml("pipeline.yml")
gtf <- ini$gtf
genome <- ini$genome
bam <- ini$bamdir

bambuAnnotations <- prepareAnnotations(gtf)
bam_dir <- bam
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Define your BAM files
bamfiles <- gsub("^\\./", "", bam_files)

# Create BamFileList using do.call
bamFiles <- BamFileList(bamfiles)
se.multiSample <- bambu(reads = bamFiles, annotations = bambuAnnotations, genome = genome)

save.dir <- "Bambu_output"
writeBambuOutput(se.multiSample, path=save.dir)

library(ggplot2)
pdf("Bambu_output/Bambu_isoform_analysis_heatmap.pdf")
plotBambu(se.multiSample, type = "heatmap")
dev.off()

pdf("Bambu_output/Bambu_isoform_analysis_PCAplot.pdf")
plotBambu(se.multiSample, type = "pca")
dev.off()
