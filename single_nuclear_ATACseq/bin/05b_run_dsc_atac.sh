#!/bin/bash

LIBRARY=$1

cd aarao_scripts/slurm/tool_specific/demuxlet

VCF_FILE=/krummellab/data1/immunox/${PROJECT}/genotypes/glands_genotyping/topmed_filt_peaks.vcf.gz  #staging/glands_genotyping/${LIBRARY}_peak_filt.vcf.gz

CR_OUT=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/cellranger/
DMX_OUT=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet/
mkdir -p ${DMX_OUT}

BARCODES=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/cell_filter/post_amulet_barcodes_of_interest_filt.list

bash spawn.sh BAMFILE=${CR_OUT}/possorted_bam.bam \
      BARCODELIST=${BARCODES} SAMPLE=${LIBRARY} OUTDIR=${DMX_OUT} VCF=${VCF_FILE} NO_TAG_UMI PARTITION=freecycle,common,krummellab DSC_ONLY

