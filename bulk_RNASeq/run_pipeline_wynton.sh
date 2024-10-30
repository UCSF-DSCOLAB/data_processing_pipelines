#!/bin/bash           # The shell language when run outside of the job scheduler

#$ -S /bin/bash         # The shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # Job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -pe smp 1          # Request 1 CPU core
#$ -l mem_free=2G     # Job requires up to 2 GiB of RAM per slot
#$ -l scratch=50G     # Job requires up to 50 GiB of local space
#$ -l h_rt=24:00:00  # Job requires up to 168 hours of runtime (7 days)
#$ -r y               # If job crashes, it should be2 restarted
#$ -o /wynton/scratch/$USER/logs/bulk_rnaseq_nf_$JOB_ID.log

## To run:
##   qsub ./run_pipeline_wynton.sh

## Arguments:
## Pass any additional arguments to Nextflow as needed (e.g., -resume, -with-timeline, -profile test)

## Create a working directory
nf_work="/wynton/scratch/$USER/nextflow/$JOB_ID/"
mkdir -p "$nf_work"
export NXF_WORK="$nf_work"
export APPTAINERENV_TMPDIR="$nf_work"
set +o posix
## Run the pipeline
nextflow run bulk_rna_seq.nf -c config/base.config -profile hpc "${@:1}"