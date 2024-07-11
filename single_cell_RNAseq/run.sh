#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/krummellab/data1/%u/logs/scseq_nf_%j.log
#SBATCH --partition=krummellab,common

# to run:
#   sbatch ./run.sh <parameter_file>.json <step>

# Arugments:
# $1 is the parameter file: e.g. example/params.json
# $2 is the step of the pipeline you are running, must be one of "pre_qc", "post_qc", "pre_fmx_qc", "post_fmx_qc"
# then pass as many additional arguments to nextflow as you'd like (e.g. -with-timeline, -profile test)


export NXF_JAVA_HOME="/krummellab/data1/erflynn/software/java/jdk-17.0.5"
export PATH=$PATH:"/krummellab/data1/erflynn/software/nextflow/22.10.4_build_5836/"

# check on arguments
failed=false
PARAM_FILE=$1
STEP=$2 

if [ ! -f "$PARAM_FILE" ]; 
then
 echo "Parameter file does not exist"
 failed=true;
fi


VALID_STEPS=("pre_qc post_qc pre_fmx_qc post_fmx_qc")
if [[ ! " ${VALID_STEPS[*]} " =~ " ${STEP} " ]];
then
  echo "please include a valid step as the second argument -- must be one of 'pre_qc', 'post_qc', 'pre_fmx_qc', 'post_fmx_qc'"
  failed=true;
fi


if [ "$failed" = true ];
then
  exit 1
else
    unset SBATCH_PARTITION
    # create a working directory
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    mkdir -p $nf_work
    export NXF_WORK=${nf_work}

    # run the pipeline
    nextflow run pipeline_${STEP}.nf -c tool.config -params-file $PARAM_FILE "${@:3}" 

fi

