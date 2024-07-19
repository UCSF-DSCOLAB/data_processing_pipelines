#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH -e /krummellab/data1/%u/logs/run_archR-%j.err
#SBATCH -o /krummellab/data1/%u/logs/run_archR-%j.out

R_CONTAINER=/krummellab/data1/singularity_images/RSingleCell/v4/RSingleCell.sif


singularity exec -B /krummellab/data1/ -B ${TMPDIR}:/tmp/ ${R_CONTAINER} Rscript 06_run_archR_post_demux.R $1

