#!/bin/bash
#SBATCH --partition=freecycle,krummellab,common
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --time=24:00:00
#SBATCH -e /krummellab/data1/%u/logs/run_amulet-%j.err
#SBATCH -o /krummellab/data1/%u/logs/run_amulet-%j.out

LIBRARY=$1
PARENT_DIR=/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed
REF_DIR=/krummellab/data1/erflynn/ref
SOFTWARE_DIR=/krummellab/data1/erflynn/software
CODE_DIR=/krummellab/data1/DSCoLab/AUTOIPI/snatac_pipeline_v0/03a_run_amulet/
RSINGULARITY_IMG=/krummellab/data1/singularity_images/RSingleCell/v3/RSingleCell.sif 

blacklist=${REF_DIR}/blacklist/hg38-blacklist.v2.bed.gz
bamfile=${PARENT_DIR}/${LIBRARY}/cellranger/possorted_bam.bam # from CR
barcodes=${PARENT_DIR}/${LIBRARY}/cell_filter/pre_amulet_barcodes_filt.csv

autosomes=${REF_DIR}/chr_list.txt # generated
output=${PARENT_DIR}/${SAMPLE}/amulet/

script_path=${SOFTWARE_DIR}/AMULET

mkdir -p $output

source ${SOFTWARE_DIR}/miniconda/etc/profile.d/conda.sh
conda activate amulet-env
${SOFTWARE_DIR}/AMULET/AMULET.sh --forcesorted --bcidx 0 --cellidx 0 --iscellidx 9 $bamfile $barcodes $autosomes $blacklist $output $script_path

singularity exec -B ${PARENT_DIR} ${RSINGULARITY_IMG} ${CODE_DIR}/post_amulet_filter.R ${LIBRARY}
