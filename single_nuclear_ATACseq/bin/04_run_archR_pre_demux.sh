#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=4:00:00
#SBATCH -e /krummellab/data1/%u/logs/run_archR-%j.err
#SBATCH -o /krummellab/data1/%u/logs/run_archR-%j.out


LIBRARY=$1
R_CONTAINER=/krummellab/data1/singularity_images/RSingleCell/v4/RSingleCell.sif
DATA_DIR=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed
CODE_DIR=${PWD}

singularity exec -B ${CODE_DIR} -B ${DATA_DIR} -B ${TMPDIR}:/tmp/ ${R_CONTAINER} \
    Rscript ${CODE_DIR}/04_run_archR_pre_demux.R ${LIBRARY}

