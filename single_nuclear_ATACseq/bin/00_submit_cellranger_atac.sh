#!/bin/bash

DATA_DIR=/krummellab/data1/immunox/${PROJECT} # CHANGE to your project
SCRIPTS_DIR=${PWD} # fill this in with the directory you copied the ATAC scripts to
#FNAMES="" # could provide this instead and use a loop

sample=$1

cd ${SCRIPTS_DIR}
export_vars="\
ALL,\
REFERENCE=/krummellab/data1/ipi/data/refs/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0,\
CELLRANGERATACVERSION=2.0.0,\
MEMORY=80"


#IFS=" "; read -a FNAMES0 <<< "$FNAMES"
#for sample in "${FNAMES0[@]}";
#do
  echo $sample
  MY_LIBRARY=${DATA_DIR}/10X_library_csvs/${sample}_libraries.csv

  echo "fastqs,sample,library_type
${DATA_DIR}/data/single_nuclear_ATAC/raw/${sample},${sample},ATAC Seq" > $MY_LIBRARY

     export_vars2="\
${export_vars},\
LIBRARIES_CSV=${MY_LIBRARY},\
SAMPLE=${sample},\
OUTDIR=${DATA_DIR}/data/single_nuclear_ATAC/processed/${sample}/cellranger"


     sbatch --export=${export_vars2} \
          --partition=freecycle,krummellab,common \
          -e /krummellab/data1/${USER}/logs/cellranger_atac_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
          -o /krummellab/data1/${USER}/logs/cellranger_atac_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
          -J cellranger_atac_count_${sample} \
          --cpus-per-task=15 \
          --mem-per-cpu=5gb \
          --time=1-00:00:00 \
          00_run_cr-atac2.0.sh
#done





