#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --time=24:00:00
#SBATCH --output=/c4/home/%u/nextflow_run_%j.log
#SBATCH --partition=krummellab,common

# Arugment:
# $1 is the parameter file: e.g. example/params.yml
# $2 is the tools file: e.g. tool.config - TODO make default
# TODO: add an option to resume

export PATH=$PATH:/krummellab/data1/ipi/software/nextflow/21.04.3_build_5560/
nextflow run pipeline_no_ribo.nf -c $2 -params-file $1 


