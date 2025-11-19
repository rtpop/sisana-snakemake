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

# Containers
SISANA_CONTAINER = config["sisana_container"]


## ----------------- ##
## SiSaNA parameters ##
## ----------------- ##

SISANA_CONFIG = config["sisana_config"]
OUTPUT_FORMAT = config["output_format"]

## Input files
MOTIF_PRIOR = os.path.join(DATA_DIR, config["motif_prior_file"])
PPI_PRIOR = os.path.join(DATA_DIR, config["ppi_prior_file"])

## output files
EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "preprocess", str(os.path.basename(EXP_SUBSET).replace(".csv", "") + "_preprocessed.txt"))
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
        LIONESS_NET, \
        LIONESS_INDEGREE, \
        LIONESS_OUTDEGREE, \
        LIONESS_EDGES, \
        if OUTPUT_FORMAT == "RData" 
            LIONESS_R
        elif OUTPUT_FORMAT == "text" 
            LIONESS_TXT