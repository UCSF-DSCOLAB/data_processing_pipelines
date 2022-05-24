#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --time=48:00:00
#SBATCH --output=/krummellab/data1/%u/logs/nextflow_run_%j.log
#SBATCH --partition=krummellab,common

# to run:
#   sbatch ./run.sh <parameter_file>.yml 

# Arugments:
# $1 is the parameter file: e.g. example/params.yml
# then pass as many additional arguments to nextflow as you'd like (e.g. -resume, -with-timeline)

function cleanup()
{
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    rm -rf ${nf_work}
}
trap cleanup EXIT


export PATH=$PATH:/krummellab/data1/ipi/software/nextflow/21.04.3_build_5560/

# check on arguments
failed=false
PARAM_FILE=$1

if [ ! -f "$PARAM_FILE" ]; 
then
 echo "Parameter file does not exist"
 failed=true;
fi


if [ "$failed" = true ];
then
  exit 1
else
    # create a working directory
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    mkdir -p $nf_work
    export NXF_WORK=${nf_work}

    # run the pipeline
    nextflow run pipeline.nf -c tool.config -params-file $PARAM_FILE "${@:2}" 
fi

