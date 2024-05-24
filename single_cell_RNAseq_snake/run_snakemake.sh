#!/bin/bash
# must have performed step 1 of USAGE in the README before running this script
module load CBI miniconda3/23.3.1-0-py39
# Load  your locally installed mamba environment for running snakemake 
conda activate depool_snakemake
python run_pipeline.py depooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/" --local-cores 2

