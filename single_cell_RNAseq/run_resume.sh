#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/krummellab/data1/%u/logs/scseq_nf_resume_%j.log
#SBATCH --partition=krummellab,common

# to run:
#   sbatch ./run_resume.sh <parameter_file>.json <step> <JOB_ID>

# Arugments:
# $1 is the parameter file: e.g. example/params.json
# $2 is the step of the pipeline you are running, must be one of "pre_qc", "post_qc", "pre_fmx_qc", "post_fmx_qc"
# $3 is the JOB_ID of the main run that you want to resume
# then pass as many additional arguments to nextflow as you'd like (e.g. -with-timeline, -profile test, etc.)


export NXF_JAVA_HOME="/krummellab/data1/erflynn/software/java/jdk-17.0.5"
export PATH=$PATH:/krummellab/data1/erflynn/software/nextflow/22.10.4_build_5836/

# check on arguments
failed=false
PARAM_FILE=$1
STEP=$2
PREV_JOB_ID=$3

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

nf_work=/c4/scratch/${USER}/nextflow/${PREV_JOB_ID}/
if [[ ! -d ${nf_work} ]];
then
echo "no previous job found to resume with job id ${PREV_JOB_ID}. please start a new run or fix the id"
failed=true;
fi


if [ "$failed" = true ];
then
  exit 1
else
    unset SBATCH_PARTITION
    # create a working directory
#    mkdir -p $nf_work
    export NXF_WORK=${nf_work}

    # run the pipeline
    nextflow run pipeline_${STEP}.nf -c tool.config -params-file $PARAM_FILE -resume "${@:4}" 
fi

