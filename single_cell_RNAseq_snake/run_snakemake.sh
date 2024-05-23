#!/bin/bash
module load mamba
source /krummellab/data1/software/miniforge3/bin/activate snakehiss
python run_pipeline.py depooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/" --local-cores 2

