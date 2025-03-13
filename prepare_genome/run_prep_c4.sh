#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/krummellab/data1/%u/logs/bulk_rnaseq_nf_%j.log
#SBATCH --exclude=c4-n20

# to run:
#   sbatch ./run_pipeline.sh -profile hpc

# Arugments:
# pass as many additional arguments to nextflow as you'd like (e.g. -resume, -with-timeline, -profile test)

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



    # create a working directory
    nf_work=/c4/scratch/${USER}/nextflow/${SLURM_JOB_ID}/
    mkdir -p $nf_work
    export NXF_WORK=${nf_work}
    export APPTAINERENV_TMPDIR=${nf_work}
    # run the pipeline

 nextflow run prepare_reference_genome.nf -c config/base.config -profile hpc_c4

