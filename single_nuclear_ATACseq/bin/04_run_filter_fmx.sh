#!/bin/bash

NUM_SAMPLES=18

# /krummellab/data1/singularity_images/bedtools/2.29.2/bedtools.sif 
# create a merged peaks.bed
MERGE_DIR="/krummellab/data1/erflynn/aip1_scratch/data/freemuxlet_data/merged_demux_AIP1-POOL-KA11-SNA"
cd ${MERGE_DIR}
DATA_DIR="/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed/"
I="AIP1-POOL-KA11-SNA1 AIP1-POOL-KA11-SNA2 AIP1-POOL-KA11-SNA3 AIP1-POOL-KA11-SNA4"
touch peaks.all.bed

IFS=" "; read -a FNAMES0 <<< "$I"
for f in "${FNAMES0[@]}"; do
    bed_file=${DATA_DIR}/${f}/cellranger/atac_peaks.bed # different for multiome
    bedtools sort -i ${bed_file} | bedtools merge -i - | cut -f 1,2,3 >> peaks.all.bed;
done

# filter the bravo reference
BRAVO_REF=/krummellab/data1/erflynn/aip1_scratch/data/dsc_filtered_ref/vcfs/sorted_bravo_dbsnp_maf0.001_pass_vrt1.vcf.gz
bcftools view -R peaks.all.bed -O z -o bravo_peaks.vcf.gz ${BRAVO_REF}

# run freemuxlet

cd /krummellab/data1/erflynn/aarao_scripts/slurm/tool_specific/freemuxlet

VCF_FILE=${MERGE_DIR}/bravo_peaks.vcf.gz
IFS=" "; read -a FNAMES0 <<< "$I"
for LIBRARY in "${FNAMES0[@]}"; do
  echo ${LIBRARY}
  CR_OUT=/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed/${LIBRARY}/cellranger/
  FMX_OUT=/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet/
  mkdir -p ${FMX_OUT}
  
  BARCODES=/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed/${LIBRARY}/cell_filter/post_amulet_barcodes_of_interest_filt.list
  
  # bam file different for multiome
  bash spawn.sh BAMFILE=${CR_OUT}/atac_possorted_bam.bam BARCODELIST=${BARCODES} SAMPLE=${LIBRARY} OUTDIR=${FMX_OUT} ONEKGENOMESVCF=${VCF_FILE} NUMSAMPLES=${NUM_SAMPLES} NO_TAG_UMI PARTITION=freecycle
done


# merge freemuxlet


# run freemuxlet
