#!/bin/bash

# load C4s conda environment
module load CBI miniconda3/23.3.1-0-py39

# create the a snakemake conda environment
conda create -n depool_snakemake python=3.11.0

#activate the environment and install mamba
conda activate depool_snakemake
conda install -c conda-forge mamba=1.5.8

# install snakemake using mamba
mamba install snakemake=8.10.7 snakemake-executor-plugin-slurm=0.4.4 -c conda-forge -c bioconda
