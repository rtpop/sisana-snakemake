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
NCORES = config["ncores"]

# Containers
SISANA_CONTAINER = config["sisana_container"]

# Directories
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
SRC = config["src_dir"]
SISANA_DIR = os.path.join(OUTPUT_DIR, "sisana")

## ----------------- ##
## SiSaNA parameters ##
## ----------------- ##

SISANA_CONFIG = config["sisana_config"]
OUTPUT_FORMAT = config["output_format"]

## Input files
MOTIF_PRIOR = os.path.join(DATA_DIR, config["motif_prior_file"])
PPI_PRIOR = os.path.join(DATA_DIR, config["ppi_prior_file"])
EXP_FILE = os.path.join(DATA_DIR, config["exp_file"])
SAMPLE_ANNO = os.path.join(DATA_DIR, config["sample_anno"])
GENE_ANNO = os.path.join(DATA_DIR, config["gene_anno"])

## output files
EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "preprocess", str(os.path.basename(EXP_FILE).replace(".csv", "") + "_preprocessed.txt"))
PANDA_NET = os.path.join(SISANA_DIR, "network", "panda_network.txt")
LIONESS_NET = os.path.join(SISANA_DIR, "network", "lioness_networks.npy")
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
        ncores = NCORES
    container:
        SISANA_CONTAINER
    message:
        "; Generating SiSaNA config file with script {params.script}"
    shell:
        """
        python {params.script} --exp {input.exp} --motif {params.motif} --ppi {params.ppi} --number {params.number} --outdir {params.outdir} --method {params.method} --pandafilepath {params.pandafilepath} --lionessfilepath {params.lionessfilepath} --output {output.config_file}
        """

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
        script = os.path.join(SRC, "utils/convert_npy_to_txt.py")
    container:
        SISANA_CONTAINER
    message:
        "; Converting LIONESS networks from npy to txt."
    shell:
        """
        python {params.script} --input {input.lioness} --output {output.lioness_txt} --pickle {input.lioness_pickle} --edges {output.lioness_edges}
        """