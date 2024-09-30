"""===========================
Pipeline dtu
By Chen-Yi Wang (2024/Sep/10)
==============================

Overview
========

This code is a pipeline for differential transcript usage (dtu) analysis, 
using a count matrix generated from Bambu


Pipeline tasks
==============

The pipeline consists of the following steps:
* filtering and dtu analysis by DRIMSeq
* statistic analysis by StageR
* final data visualisation


Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin dtu config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin dtu make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin dtu make full -v5 --local

Input files
===========

Input files should be a count matrix of transcripts output by Bambu

Pipeline output
===============

The pipeline outputs include count matrices of results and significant genes and the associated transcripts,
and plots for visualization


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
RMD_ROOT = os.path.join(os.path.dirname(__file__), "pipeline_dtu")
# R folder in main directory

R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

@follows(mkdir("DTU_DRIMSeq.dir"))
def dtu():
    '''
    R script task to run DRIMSeq and StageR
    '''
    R_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"R"))

    job_memory = "70G"

    statement = '''
    Rscript %(R_PATH)s/DTU_DRIMSeq.R '''

    P.run(statement)


@follows(dtu)
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
