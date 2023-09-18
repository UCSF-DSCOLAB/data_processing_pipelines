#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=72:00:00
#SBATCH --output=/krummellab/data1/%u/logs/bulk_rnaseq_nf_%j.log
#SBATCH --partition=krummellab,common

# to run:
#   sbatch ./run.sh <parameter_file>.json <step>

# Arugments:
# pass as many additional arguments to nextflow as you'd like (e.g. -with-timeline, -profile test)

function cleanup()
{
    if [ $? -eq 0 ]
    then
    echo "step complete done, deleting working directory"
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    rm -rf ${nf_work}
    fi
}
trap cleanup EXIT


# check on arguments
failed=false


if [ "$failed" = true ];
then
  exit 1
else
    # create a working directory
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    mkdir -p $nf_work
    export NXF_WORK=${nf_work}
    export APPTAINERENV_TMPDIR=${nf_work}
    # run the pipeline
    nextflow run bulk_rna_seq.nf "${@:1}"
fi
