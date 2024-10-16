# **Bulk ONT Transcriptomics Analysis Workflow**

This repository provides workflows for bulk Oxford Nanopore Technologies (ONT) transcriptomics analysis, designed to streamline the entire process using `fastq.gz` files as input. The input data is processed through the following steps:

1. **Pychopper**: For trimming and identifying full-length transcripts.
2. **Minimap2**: For mapping reads to the reference genome.
3. **Bambu**: To generate gene- or transcript-level count matrices.

The pipeline also includes downstream analyses, such as:

  - Differential transcript usage (DTU) analysis
  - Protein domain search
  - Protein alignment

# **Installation**

The workflows are structured using CGAT's Ruffus pipeline framework and are currently built within the TallyTrin environment. For detailed installation instructions, please refer to [TallyTrin](https://github.com/cribbslab/TallyTriN).

For the R packages used in the workflows, please visit their respective GitHub pages for installation and further details.

# **Requirements and Reference Genome**

Please check the `pipeline.yml` file of each pipeline for information on software requirements and the genome references needed for data processing. 
