#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=4:00:00
#SBATCH -e /krummellab/data1/%u/logs/run_archr-%j.err
#SBATCH -o /krummellab/data1/%u/logs/run_archr-%j.out

CODE_DIR=$PWD
singularity exec -B /krummellab/data1/immunox/${PROJECT}/ \
                 -B /krummellab/data1/DSCoLab/${PROJECT}/ \
                 -B ${TMPDIR}:/tmp/ \
                 /krummellab/data1/singularity_images/RStudioSingleCell/v4/RStudioSingleCell.sif Rscript \
    ${CODE_DIR}/01_load.R ${1}
