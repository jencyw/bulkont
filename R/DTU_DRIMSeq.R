
### Author: Chen-Yi Wang
### Title: "Differential transcript usage (DTU) analysis using DRIMSeq and StageR"
### Goal: to perform DTU analysis and visualize DTU at the detected gene
### Output: results of DTU at transcript and gene level, 
### Output: data visualization includes lineplot, ribbon plot, boxplot, heatmap and transcript structures on conditions
### Before run, please check the following requirement:   
# required files: sample_sheet.csv, count matrix, and annotation reference
# R availability in the environment
# Availability of the required R packages

library(data.table)
library(DEXSeq)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(dplyr)
library(DRIMSeq)
library(pheatmap)
library(yaml)

#devtools::install_github("dzhang32/ggtranscript")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("stageR")

ini <- read_yaml("pipeline.yml")
min_samps_feature_expr <- ini$min_samps_feature_expr
min_feature_expr <- ini$min_feature_expr
min_samps_gene_expr <- ini$min_samps_gene_expr
min_gene_expr <- ini$min_gene_expr
padj <- ini$padj
control <- ini$control
treatment <- ini$treatment

### Load data and annotation reference
# sample sheet: sample_id and condition
samps <- read.table("sample_sheet.csv", sep = ",", header = TRUE)
samps <- samps %>%
  filter(condition==control|condition==treatment)
samps = data.frame(sample_id = samps$sample_id, group = samps$condition)

# Load count matrix
# make sure the order of the columns are the same with the order in the sample sheet
df <- read.table("Bambu_output/counts_transcript.txt", sep = "\t", header = TRUE)
colnames(df) <- gsub("\\.bam$", "", colnames(df)) 
df <- df[!grepl("^Bambu", df$TXNAME),]
txdf <- df[,c(1:2)]

rownames(df)<- df$TXNAME

cts <- df %>%
  dplyr::select(any_of(samps$sample_id))
cts <- cts[rowSums(cts)> 0,]

txdf$TXNAME_2 <- gsub("\\..*", "",txdf$TXNAME) 
txdf$ensembl_gene_id <- gsub("\\..*", "",txdf$GENEID)

txdf.sub = txdf[match(rownames(cts),txdf$TXNAME),]
counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME, cts)

# DRIMSeq analysis
d = DRIMSeq::dmDSdata(counts = counts, samples = samps)
d <- DRIMSeq::dmFilter(d,
                       min_samps_feature_expr=min_samps_feature_expr, min_feature_expr=min_feature_expr,
                       min_samps_gene_expr=min_samps_gene_expr, min_gene_expr=min_gene_expr)
design_full <- model.matrix(~ group, data = samples(d))

dir.create("DTU_DRIMSeq.dir")
set.seed(123)
## Calculate precision
d <- dmPrecision(d, design = design_full)
common_precision(d)
plotPrecision(d)
ggp <- plotPrecision(d)
pdf("DTU_DRIMSeq.dir/plotPrecision.pdf")
ggp + geom_point(size = 4)
dev.off()

coef <- paste0("group", treatment)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = coef, verbose = 1)

design(d)

pdf("DTU_DRIMSeq.dir/plotPValues_feature.pdf")
plotPValues(d, level = "feature")
dev.off()

strp <- function(x) substr(x,1,15)
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]
res$ensembl_gene_id <- strp(res$gene_id)
no.na <- function(x) ifelse(is.na(x), 1, x)
library(biomaRt)
mart <- useMart("ensembl","hsapiens_gene_ensembl")##hsapiens_gene_ensembl
ensemble2gene <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = res$ensembl_gene_id, 
                       mart = mart)
res_gene <- merge(res, ensemble2gene,by="ensembl_gene_id", all.x=T)
res_gene$pvalue <- no.na(res_gene$pvalue)
res_gene$adj_pvalue <- no.na(res_gene$adj_pvalue)
write_tsv(res_gene, file = "DTU_DRIMSeq.dir/res_genes.tsv")

# p-value results at the transcript level: to test whether the proportions of the transcript changed within the gene
res.txp <- DRIMSeq::results(d, level="feature")
res.txp$pvalue <- no.na(res.txp$pvalue)
res.txp$adj_pvalue <- no.na(res.txp$adj_pvalue)
res.txp$lr <- no.na(res.txp$lr)
res.txp$ensembl_transcript_id <- strp(res.txp$feature_id)

txp.info <- read_tsv(txp_info)
res.txp.info <- merge(res.txp, txp.info, by.x= "ensembl_transcript_id", by.y="tx_id", all.x=T)
write_tsv(res.txp.info, file = "DTU_DRIMSeq.dir/res_transcripts_info_Ensemblv110.tsv")

idx <- which(res$adj_pvalue < padj)
sig <- res[idx,]
sig_symbol <- merge(sig, ensemble2gene, by ="ensembl_gene_id")
sig_symbol <- sig_symbol[!duplicated(sig_symbol$ensembl_gene_id),]
sig_symbol <- sig_symbol%>%
  filter(sig_symbol$external_gene_name!="")
filename <- paste0("DTU_DRIMSeq.dir/DRIMSeq_sig_padj", padj, "_biomart.tsv")
write_tsv(sig_symbol, file = filename)

### StageR 
library(stageR)
smallProportionSD <- function(d, filter = 0.1) {
  # Generate count table
  cts = as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
  # Summarise count total per gene
  gene.cts = rowsum(cts, counts(d)$gene_id)
  # Use total count per gene as count per transcript
  total.cts = gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  # Calculate proportion of transcript in gene
  props = cts/total.cts
  rownames(props) = rownames(total.cts)
  
  # Calculate standard deviation
  propSD = sqrt(rowVars(props))
  # Check if standard deviation of per-sample proportions is < 0.1
  propSD < filter
}

filt = smallProportionSD(d)

res.txp.filt = DRIMSeq::results(d, level = "feature")
res.txp.filt$pvalue[filt] = 1
res.txp.filt$adj_pvalue[filt] = 1

pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)
pConfirmation <- matrix(res.txp.filt$pvalue, ncol=1)
dimnames(pConfirmation) <- list(strp(res.txp.filt$feature_id), "transcript")
tx2gene <- data.frame(res.txp[,c("feature_id", "gene_id")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene[,1:2])
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})
drim_padj_symbol <- merge(drim.padj, txp.info, by.x= "txID", by.y="tx_id", all.x=T)
write_tsv(drim_padj_symbol, file = "DTU_DRIMSeq.dir/StageR_sig_padj0.05_Ensemblv110.tsv")

### Lineplot of genes with DTU

dir.create("DTU_DRIMSeq.dir/StageR_lineplot_DTU")
df_plot <- drim_padj_symbol%>%
  filter(!duplicated(geneID))%>%
  filter(external_gene_name!="NA")
df_plot <- merge(df_plot, res, by.x="geneID", by.y="ensembl_gene_id")
gene_names <- df_plot$gene_id
names(gene_names) <- df_plot$external_gene_name

for (alias in names(gene_names)){
  
  gene_name <- gene_names[alias]
  #Generate plot for the current gene
  # Replace `your_plot_function` with the function you're using to generate the plot
  plot <- plotProportions(d, gene_id = gene_name, group_variable = "group", group_colors= c("#293f6e", "#f96866"), plot_type = "lineplot")
  
  # Construct filename for the plot
  # You may need to adjust the directory path as per your requirement
  filename <- paste0("DTU_DRIMSeq.dir/StageR_lineplot_DTU/", alias, ".pdf")
  
  # Save the plot
  ggsave(filename, plot, width = 10, height = 6)
}

### Ribbon plots of genes with DTU
dir.create("DTU_DRIMSeq.dir/StageR_ribbonplot_DTU")

for (alias in names(gene_names)){
  
  gene_name <- gene_names[alias]
  
  # Generate plot for the current gene
  # plot
  plot <- plotProportions(d, gene_id = gene_name, group_variable = "group",plot_type = "ribbonplot")
  
  filename <- paste0("DTU_DRIMSeq.dir/StageR_ribbonplot_DTU/", alias, ".pdf")
  
  # Save the plot
  ggsave(filename, plot, width = 8, height = 6)
}

### boxplot of genes with DTU

dir.create("DTU_DRIMSeq.dir/StageR_boxplot_condition")

for (alias in names(gene_names)){
  
  gene_name <- gene_names[alias]
  #Generate plot for the current gene
  # Replace `your_plot_function` with the function you're using to generate the plot
  plot <- plotProportions(d, gene_id = gene_name, group_variable = "group", group_colors= c("#293f6e", "#f96866"), plot_type = "boxplot1")
  
  filename <- paste0("DTU_DRIMSeq.dir/StageR_boxplot_condition/", alias, ".pdf")
  
  # Save the plot
  ggsave(filename, plot, width = 10, height = 6)
}

dir.create("DTU_DRIMSeq.dir/StageR_boxplot_transcript")

for (alias in names(gene_names)){
  
  gene_name <- gene_names[alias]
  #Generate plot for the current gene
  # Replace `your_plot_function` with the function you're using to generate the plot
  plot <- plotProportions(d, gene_id = gene_name, group_variable = "group",plot_type = "boxplot2")
  
  filename <- paste0("DTU_DRIMSeq.dir/StageR_boxplot_transcript/", alias, ".pdf")
  
  # Save the plot
  ggsave(filename, plot, width = 14, height = 6)
}

# Heatmap of gene with DTU
dir.create("DTU_DRIMSeq.dir/StageR_df_plotproportions")
dir.create("DTU_DRIMSeq.dir/StageR_Heatmap")

for (alias in names(gene_names)){
  
  gene_name <- gene_names[alias]
  #Generate plot for the current gene
  # Replace `your_plot_function` with the function you're using to generate the plot
  plot <- plotProportions(d, gene_id = gene_name, group_variable = "group", group_colors= c("#293f6e", "#f96866"), plot_type = "lineplot")
  plotmatrix <- as.data.frame(plot$plot_env$prop_samp)
  plotmatrix$ensembl_transcript_id <- gsub("\\..*", "",plotmatrix$feature_id)
  plotmatrix2 <- merge(plotmatrix, txp.info, by.x= "ensembl_transcript_id", by.y="tx_id", all.x=T)
  plotmatrix2$tx_character <- paste0(plotmatrix2$ensembl_transcript_id, "_", plotmatrix2$tx_biotype, "_", plotmatrix2$protein_aa_length, "aa")
  
  # Reshape the data for heatmap
  data_wide <- plotmatrix2[,c(1:5,29)]
  data_wide <- dcast(data_wide, tx_character ~ sample_id, value.var="proportion",fun.aggregate = sum)
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(data_wide[, -1])
  rownames(heatmap_matrix) <- data_wide$tx_character
  
  # Plot the heatmap
  heatmap <- pheatmap(heatmap_matrix, cluster_rows = TRUE, cluster_cols = FALSE, display_numbers = TRUE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
  
  
  # Construct filename for the plot
  # You may need to adjust the directory path as per your requirement
  filename1 <- paste0("DTU_DRIMSeq.dir/StageR_df_plotproportions/df_",alias, ".tsv")
  filename2 <- paste0("DTU_DRIMSeq.dir/StageR_Heatmap/", alias, ".pdf")
  
  # Save the datamatrix
  write_tsv(plotmatrix2, filename1)
  # Save the plot
  ggsave(filename2, heatmap, width=10, height = 8)
}

# visualization of transcript structure and annotation
dir.create("DTU_DRIMSeq.dir/StageR_transcript_structure_annotation")
library(magrittr)
library(ggtranscript)
library(rtracklayer)
gtf <- read_tsv(gtf)

for (alias in names(gene_names)){
  
  gene_of_interest <- alias
  
  annotation_from_gtf <- gtf %>% 
    dplyr::filter(
      !is.na(gene_name), 
      gene_name == gene_of_interest
    ) 
  
  # extract the required annotation columns
  annotation_from_gtf <- annotation_from_gtf %>% 
    dplyr::select(
      seqid,
      start,
      end,
      strand,
      type,
      gene_name,
      transcript_id,
      transcript_name,
      transcript_biotype
    )
  
  # extract exons
  exons <- annotation_from_gtf %>% dplyr::filter(type == "exon")
  
  plot <- exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(
      aes(fill = transcript_biotype)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id"),
      aes(strand = strand), 
      arrow.min.intron.length = 2500
    )
  
  # Construct filename for the plot
  # You may need to adjust the directory path as per your requirement
  filename <- paste0("DTU_DRIMSeq.dir/StageR_transcript_structure_annotation/", alias, ".pdf")
  
  # Save the plot
  ggsave(filename, plot, width = 10, height = 6)
}

### Pathway Analysis of significant genes in DRIMSeq
dir.create("DTU_DRIMSeq.dir/Pathway_analysis")
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
library <- as.data.frame(dbs$libraryName)

dbs <- c("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2021", "GO_Biological_Process_2023", "MSigDB_Hallmark_2020", 	
         "MSigDB_Oncogenic_Signatures","Panther_2016", "WikiPathway_2023_Human")
genelist <- drim_padj_symbol$external_gene_name
websiteLive <- TRUE

if (websiteLive) {
  enriched <- enrichr(genelist, dbs)
}

MsigDB_Hallmark <- as.data.frame(subset(enriched[["MSigDB_Hallmark_2020"]], P.value < 0.05))
KEGG <- as.data.frame(subset(enriched[["KEGG_2021_Human"]], P.value < 0.05))
wiki <- as.data.frame(subset(enriched[["WikiPathway_2023_Human"]], P.value < 0.05))
Reactome <- as.data.frame(subset(enriched[["Reactome_2022"]], P.value < 0.05))
GOBP_2021 <- as.data.frame(subset(enriched[["GO_Biological_Process_2021"]], P.value < 0.05))
GOBP_2023 <- as.data.frame(subset(enriched[["GO_Biological_Process_2023"]], P.value < 0.05))
Panther <- as.data.frame(subset(enriched[["Panther_2016"]], P.value < 0.05))

write_tsv(MsigDB_Hallmark, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_MsigDB_Hallmark.tsv")
write_tsv(KEGG, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_KEGG.tsv")
write_tsv(wiki, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_wiki.tsv")
write_tsv(Reactome, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_Reactome.tsv")
write_tsv(GOBP_2021, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_GOBP_2021.tsv")
write_tsv(GOBP_2023, file = "DTU_DRIMSeq.dir/Pathway_analysis/StageR_GOBP_2023.tsv")

matrix_names  <- c("MsigDB_Hallmark", "KEGG","wiki","Reactome","GOBP_2021","GOBP_2023")

for (matrix_name in matrix_names){
  
  matrix <- get(matrix_name)
  matrix$rank<- -log10(as.numeric(matrix$P.value))
  positions <- matrix$Term [1:20]
  
  plot <- ggplot(matrix,aes(Term, rank))+geom_bar(stat = "identity",fill="#2f7ce8",width = 0.6)+
    theme_bw() + 
    geom_blank()+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  plot <- plot+
    scale_x_discrete(limits = positions)+
    theme(axis.text.y = element_text(colour = "black", size = 10))+
    coord_flip()+ 
    labs(x= "Top Enriched Pathways", y = "-log(p-value)", title = paste0("Enriched pathways in ", matrix_name))+
    scale_x_discrete(limit=rev(positions),labels = function(x) lapply(strwrap(x, width = 50, simplify = FALSE), paste, collapse="\n"))
  
  filename <- paste0("DTU_DRIMSeq.dir/Pathway_analysis/", "StageR_",matrix_name, ".pdf")
  
  # Save the plot
  ggsave(filename, plot)
}

