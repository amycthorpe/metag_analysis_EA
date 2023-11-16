##############################
# MODULES
import os, re
import glob
import pandas as pd

# ##############################
# # Parameters
# CORES=int(os.environ.get("CORES", 4))


##############################
# Paths
SRC_DIR = srcdir("../scripts")
ENV_DIR = srcdir("../envs")
NOTES_DIR = srcdir("../notes")
SUBMODULES= srcdir("../../submodules")

##############################
# Dependencies 


##############################
# default executable for snakemake
shell.executable("bash")


##############################
# working directory
workdir:
    config["work_dir"]


##############################
# Relevant directories
DATA_DIR = config["data_dir"]
DB_DIR = config["kraken2"]["db"]
RESULTS_DIR = config["results_dir"]
ENV_DIR=config['env_dir']
SRC_DIR=config['src_dir']


##############################
# Steps
STEPS = config["steps"]


##############################
# Input
# SAMPLES = [line.strip() for line in open("config/socd_sample_list").readlines()]
# Sample table (tab-separated, w/ header, 1st column is sample ID)
SAMPLES = pd.read_csv(config["samples"], header=0, sep="\t").set_index("Sample_ID", drop=False)
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]  
