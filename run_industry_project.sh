#!/bin/bash

# Source the Conda environment
source /work/cmcc/ad05922/conda_startup
conda activate pypsa-eur

# Submit the job to the batch system
bsub -q p_long -M 32G -n 1 -P 0607 snakemake -call all --configfile=config/industry_config.yaml
