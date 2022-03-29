#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --time=24:00:00
#SBATCH --output=/krummellab/data1/%u/logs/nextflow_run0_%j.log
#SBATCH --partition=krummellab,common

# Arugment:
# $1 is the parameter file: e.g. example/params.yml


export PATH=$PATH:/krummellab/data1/ipi/software/nextflow/21.04.3_build_5560/

nxf_work=/c4/scratch/nf_tmp/${USER}/${SLURM_JOB_ID}/ 
mkdir -p $nxf_work
export NXF_WORK=${nxf_work}

nextflow run pipeline.nf -params-file $1 

