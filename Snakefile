## ------------------------------------------------------------------------------------------- ##
## HOW TO RUN                                                                                  ##
## Run from directory containing Snakefile                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For dry run                                                                                 ##
## snakemake --cores 1 -np                                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For local run                                                                               ##
## snakemake --cores 1                                                                         ##
## ------------------------------------------------------------------------------------------- ##
## For running with singularity container                                                      ##
## snakemake --cores 1 --use-singularity --singularity-args '\-e' --cores 1                    ##
## ------------------------------------------------------------------------------------------- ## ------------------------------ ##
## For running in the background                                                                                                 ##
## nohup snakemake --cores 1 --use-singularity --singularity-args '\-e' 2>&1 > logs/snakemake_$(date +'%Y-%m-%d_%H-%M-%S').log & ##
## ----------------------------------------------------------------------------------------------------------------------------- ##
## How to manually run a container on biotin4 cluster with binding of the Snakemake folder but no software env                   ##
## apptainer shell --cleanenv --containall --no-home --bind path/to/snakemake/folder:/home/ .snakemake/singularity/your_image    ##
## ----------------------------------------------------------------------------------------------------------------------------- ##

##-----------##
## Libraries ##
##-----------##

import os 
import sys
import glob
from pathlib import Path
import time

## ----------------- ##
## Global parameters ##
## ----------------- ##

# Config file
global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH
NCORES = config.get("ncores", 1)

# Containers
SISANA_CONTAINER = config.get("sisana_container", "docker://rtpop/sisana_container:latest")

# Directories
DATA_DIR = config.get("data_dir", "data")
OUTPUT_DIR = config.get("output_dir", "output")
SRC = config.get("src_dir", "src")
SISANA_DIR = os.path.join(OUTPUT_DIR, "sisana")

## ----------------- ##
## SiSaNA parameters ##
## ----------------- ##

SISANA_CONFIG = config.get("sisana_config", "sisana_config.yaml")
OUTPUT_FORMAT = config.get("output_format", "RData")
INPUT_FORMAT = config.get("input_format", "txt")

## Input files
MOTIF_PRIOR = os.path.join(DATA_DIR, config.get("motif_prior_file", "motif_prior.txt"))
PPI_PRIOR = os.path.join(DATA_DIR, config.get("ppi_prior_file", "ppi_prior.txt"))
EXP_FILE = os.path.join(DATA_DIR, config.get("exp_file", "exp.txt"))
SAMPLE_ANNO = os.path.join(DATA_DIR, config.get("sample_anno", "sample_anno.txt"))
GENE_ANNO = os.path.join(DATA_DIR, config.get("gene_anno", "gene_anno.txt"))
## output files
EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "preprocess", str(os.path.basename(EXP_FILE).replace(".txt", "") + "_preprocessed.txt"))
PANDA_NET = os.path.join(SISANA_DIR, "network", "panda_network.txt")
LIONESS_FILE = config.get("lioness_file", "lioness_network.npy")
LIONESS_NET = os.path.join(SISANA_DIR, "network", LIONESS_FILE)
LIONESS_SAMPLES = os.path.join("tmp/samples.txt") # since sisana currently stores it here instead of the main output folder. Hopefully will be fixed with next update
LIONESS_OUTDEGREE = os.path.join(SISANA_DIR, "network", "lioness_outdegree.csv")
LIONESS_INDEGREE = os.path.join(SISANA_DIR, "network", "lioness_indegree.csv")

## Intermediate files
LIONESS_R = LIONESS_NET.replace(".npy", ".RData") # convert lioness results to R format
LIONESS_TXT = LIONESS_NET.replace(".npy", ".txt") # intermediate step
LIONESS_EDGES = LIONESS_NET.replace(".npy", "_edges.txt") # intermediate step
LIONESS_PICKLE = "tmp/lioness.pickle" # intermediate step

## ------ ##
## Rules  ##
## ------ ##

## Rule ALL ##
rule all:
    input:
        [
            LIONESS_NET,
            LIONESS_INDEGREE,
            LIONESS_OUTDEGREE,
            LIONESS_EDGES
        ] + ([LIONESS_R] if OUTPUT_FORMAT == "RData" else [LIONESS_TXT, LIONESS_PICKLE] if OUTPUT_FORMAT == "text" else [])

## --------------- ##
## Making networks ##
## --------------- ##

rule generate_sisana_params:
    """
    This rule generates the SiSaNA config file.
    """
    input:
        exp = EXP_FILE
    output:
        config_file = SISANA_CONFIG
    params:
        script = os.path.join(SRC, "generate_sisana_config.py"), \
        motif = MOTIF_PRIOR, \
        ppi = PPI_PRIOR, \
        number = config["min_exp_samples"], \
        outdir = os.path.join(SISANA_DIR, "preprocess"), \
        method = config["sisana_method"], \
        pandafilepath = PANDA_NET, \
        lionessfilepath = LIONESS_NET, \
        ncores = NCORES, \
        input_format = INPUT_FORMAT
    container:
        SISANA_CONTAINER
    message:
        "; Generating SiSaNA config file with script {params.script}"
    script:
        "{params.script}"

rule generate_networks:
    """
    This rule runs PANDA/LIONESS using the SiSaNA pipeline.

    SiSaNA is available at
    https://github.com/kuijjerlab/sisana

    Inputs
    ------
    SISANA_CONFIG:
        Config yml file for SiSaNA.
    ------
    Outputs
    -------
    EXPRESSION_FILTERED:
        A TXT file with filtered expression.
    MOTIF_PRIOR_FILTERED:
        A TXT file with filtered motif prior.
    PANDA_NET:
        A TXT file with the PANDA network.
    SAMPLES:
        A TXT file with the sample IDs.
    """
    input:
        sisana_config = SISANA_CONFIG
    output:
        EXPRESSION_FILTERED, \
        PANDA_NET, \
        LIONESS_NET, \
        LIONESS_SAMPLES, \
        LIONESS_PICKLE, \
        LIONESS_INDEGREE, \
        LIONESS_OUTDEGREE
    container:
        SISANA_CONTAINER
    message: 
        "; Running sisana preprocess on {input.sisana_config}."
    shell:
        """
        echo sisana preprocess {input.sisana_config}
        sisana preprocess {input.sisana_config}

        echo sisana generate {input.sisana_config}
        sisana generate {input.sisana_config}
        """

## ------------------------------------ ##
## Convert LIONESS networks to R format ##
## ------------------------------------ ##

rule convert_to_txt:
    input:
        lioness = LIONESS_NET, \
        lioness_pickle = LIONESS_PICKLE
    output:
        lioness_txt = LIONESS_TXT, \
        lioness_edges = LIONESS_EDGES
    params:
        script = os.path.join(SRC, "convert_npy_to_txt.py")
    container:
        SISANA_CONTAINER
    message:
        "; Converting LIONESS networks from npy to txt."
    script:
        "{params.script}"