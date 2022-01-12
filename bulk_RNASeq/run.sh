#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --time=24:00:00
#SBATCH --output=/krummellab/data1/%u/logs/nextflow_run_%j.log
#SBATCH --partition=krummellab,common

# Arugment:
# $1 is the parameter file: e.g. example/params.yml
# $2 is the tools file: e.g. tool.config - TODO make default
# pass as many additional arguments to nextflow as you'd like (e.g. -resume)

export PATH=$PATH:/krummellab/data1/ipi/software/nextflow/21.04.3_build_5560/
nextflow run pipeline.nf -c $2 -params-file $1 "${@:3}"

#nextflow run /krummellab/data1/erflynn/data_processing_pipelines/bulk_RNASeq/pipeline.nf -c $2 -params-file $1 -with-report $3 -with-trace -with-timeline $4


