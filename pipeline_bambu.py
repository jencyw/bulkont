"""===========================
Pipeline bambu
===========================

Overview
========

This code is a tool for transcript identification and discovery


Pipeline tasks
==============

The pipeline consists of the following steps:
* transcript identification and discovery
* output of count matrices


Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin bambu config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin bambu make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin bambu make full -v5 --local

Input files
===========

Input files should be sorted bam files of long-read sequencing

Pipeline output
===============

The pipeline outputs several counts matrices with samples as columns and
rows as either transcripts or genes.


Code
====

"""
import sys
import os
import pysam
import glob
import pandas as pd
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

# Root of Rmarkdown folder in pipeline folder
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_bambu")
# R folder in main directory

R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

@follows(mkdir("Bambu_output"))
def bambu():
    '''
    R script task to run Bambu
    '''
    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

    job_memory = "70G"

    statement = '''
    Rscript %(R_PATH)s/Bambu.R '''

    P.run(statement)


@follows(bambu)
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
