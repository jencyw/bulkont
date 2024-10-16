"""===========================
Pipeline count
===========================

Overview
========

This code is a CGAT pipeline for processing bulk ONT-transcriptomes from raw fastq file,
trimmed, mapped and performed transcript identification and discovery using Bambu. The
output count matrices at gene and transcript level are further analysed using either DESeq2 or
several downstream analyses of differential transcrip usage. The pipeline makes use of multiple
libraries and tools like pychopper, minimap2, and samtools.

Pipeline tasks
==============

The pipeline consists of the following steps:
    * Detecting UMI sequences.
    * Trimming Nanopore adaptors and UMI from reads.
    * Concatenating two kinds of fastq files output by pychopper: Full-length and rescued.
    * Map reads to transcripts using minimap2.
    * Processing and sorting BAM files using samtools.
    * Transcript identification and discovery using Bambu.
    * Diffenrential transcript usage analysis by DRIMSeq and StageR (pipeline_dtu.py).
    * Downstream analyses including protein alignment, domain search (pipeline_isoprot.py). 

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin bulkont config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin bulkont make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin bulkont make full -v5 --local

Input files
===========

Input files should be fastq.gz files of nanopore reads in the folder named "data.dir".

Pipeline output
===============

The pipeline outputs a counts matrix with sample as columns and
rows as either transcripts or genes.

Code
====

"""
import sys
import os
import pysam
import glob
import pandas as pd
import subprocess as P
from ruffus import *
import cgatcore.iotools as iotools
import cgatcore.pipeline as P
import cgatcore.experiment as E
from cgatcore.pipeline import cluster_runnable

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


SEQUENCESUFFIXES = ("*.fastq.gz")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

# R folder in main directory
R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

@follows(mkdir("pychopper_processed_fastq.dir"))
@transform(FASTQTARGET,
           regex(r"data.dir/(\S+).fastq.gz"),
           [r"pychopper_processed_fastq.dir/\1_full_length.fastq",
           r"pychopper_processed_fastq.dir/\1_rescued.fastq"])

def run_pychopper(infile, outfiles):
    '''
    This function runs pychopper on the input fastq file and producing:
    - full_length.fastq
    - rescued.fastq
    '''

    full_length_output, rescued = outfiles
    #basename = os.path.splitext(os.path.basename(infile))[0]

    statement = '''pychopper -U -y -k PCB111 -w %(rescued)s %(infile)s %(full_length_output)s'''

    P.run(statement, job_options='-t 24:00:00', log_file="pychopper.log")

@follows(mkdir("merged_fastq.dir"))
@follows(run_pychopper)
@collate(run_pychopper,
         regex(r"pychopper_processed_fastq.dir/(\S+)_full_length.fastq"),
         r"merged_fastq.dir/\1_merged.fastq")
 
def cat_fastq(infiles, outfile):
    '''
    This function concatenates the full-length and rescued fastq files
    for each individual sample.
    '''
    # Flatten the list of infiles (in case infiles is a list of lists)
    infiles_flat = [item for sublist in infiles for item in sublist]

    # Join the input files into a single string for the shell command
    infiles_str = " ".join(infiles_flat)

    # Concatenate full-length and rescued fastq files into one output file
    statement = f"cat {infiles_str} > {outfile}"

    P.run(statement, job_options='-t 24:00:00', log_file="cat_fastq.log")

@transform(cat_fastq,
           regex(r"merged_fastq.dir/(\S+).fastq"),
           r"merged_fastq.dir/\1.fastq.gz")

def gzip(infile, outfile):
    '''
    This function performs gzip.
    '''

    statement = '''gzip %(infile)s'''

    P.run(statement, job_memory="60G", job_options='-t 72:00:00',log_file="minimap2.log")


@follows(mkdir("mapped_files.dir"))
@transform(gzip,
           regex(r"merged_fastq.dir/(\S+)_merged.fastq.gz"),
           r"mapped_files.dir/\1_mapped.sam")
def run_minimap2(infile, outfile):
    '''
    This function performs mapping using minimap2.
    '''

    statement = '''minimap2 -y -ax splice --MD %(genome_fasta)s %(infile)s > %(outfile)s'''

    P.run(statement, job_memory="60G", job_options='-t 72:00:00',log_file="minimap2.log")

@transform(run_minimap2,
           regex("mapped_files.dir/(\S+)_mapped.sam"),
           r"mapped_files.dir/\1_sorted.bam")
def samtools_sort(infile, outfile):
    '''
    This function sorts the input BAM file using samtools and
    writes the output to the specified outfile.
    '''

    name = infile.replace("_mapped.sam", "")

    statement = '''samtools view -bS %(name)s_mapped.sam > %(name)s.bam &&
                   samtools sort %(name)s.bam -o %(name)s_sorted.bam &&
                   samtools index %(name)s_sorted.bam'''

    P.run(statement, job_options='-t 24:00:00',log_file="samtools_sort.log")

@follows(mkdir("Bambu_output.dir"))
@follows(samtools_sort)
def bambu():
    '''
    R script task to run Bambu
    '''
    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

    job_memory = "70G"

    statement = '''
    Rscript %(R_PATH)s/Bambu.R '''

    P.run(statement, job_options="-t 24:00:00",log_file="bambu.log")

@follows(run_pychopper, cat_fastq, gzip, run_minimap2, samtools_sort, bambu)
def full():
    '''
    A placeholder function that serves as a checkpoint
    to run all previous ruffus tasks and ensure that all
    previous tasks are completed.
    '''
    pass


def main(argv=None):
    '''
    The main function that runs the pipeline using the cgatcore.pipeline module.
    Takes an optional argument list (default is sys.argv).

    Please note that some of these functions use external Python scripts or
    tools. For a complete understanding of their functionality, it is
    necessary to examine the code of those scripts as well.
    '''
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))



