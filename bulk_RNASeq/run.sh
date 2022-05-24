#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --time=48:00:00
#SBATCH --output=/krummellab/data1/%u/logs/nextflow_run_%j.log
#SBATCH --partition=krummellab,common

# to run:
#   sbatch ./run.sh <parameter_file>.yml [keep_work|remove_work]

# Arugments:
# $1 is the parameter file: e.g. example/params.yml
# $2 is whether to remove work directory ("keep_work" or "remove_work")
# then pass as many additional arguments to nextflow as you'd like (e.g. -resume, -with-timeline)

export PATH=$PATH:/krummellab/data1/ipi/software/nextflow/21.04.3_build_5560/

# check on arguments
failed=false
NUM_ARGS=$#

if [ $NUM_ARGS -lt 2 ];
then
 echo "Please run with at least two arguments. Usage:  sbatch ./run.sh <parameter_file>.yml [keep_work|remove_work]"
 failed=true;
fi

PARAM_FILE=$1
WORK_ARG=$2

if [ ! -f "$PARAM_FILE" ]; 
then
 echo "Parameter file does not exist"
 failed=true;
fi

if [ "$WORK_ARG" != "remove_work" ] && [ "$WORK_ARG" != "keep_work" ];
then
  echo "The second argument must be either 'remove_work' or 'keep_work'"
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
    nextflow run pipeline.nf -c tool.config -params-file $PARAM_FILE "${@:3}" 


    # clean up if desired
    if [ ${WORK_ARG} == "remove_work" ]
    then
	rm -rf ${nf_work}
    fi
fi
