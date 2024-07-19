#!/bin/bash


REF_DIR=/krummellab/data1/immunox/refs/GRCh38_p13_files/bravo_dbsnp/
IMX_DIR=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/

SAMPLE_LIST=() # TODO fill in samples

# copy over all of the peak files from that tissue
for i in "${SAMPLE_LIST[@]}";
do
cp ${IMX_DIR}/${sample}/cellranger/peaks.bed ${sample}.bed
done

# merge the peak files
I=`ls *.bed`
touch peaks.all.bed

IFS=" "; read -a FNAMES0 <<< "$I"
for f in "${FNAMES0[@]}"; do
    bedtools sort -i $f | bedtools merge -i - | cut -f 1,2,3 >> peaks.all.bed;
done


# use the merged peak file to extract only the snATAC regions from the reference
bcftools view -R peaks.all.bed -O z -o topmed_filt_peaks.vcf.gz ${REF_DIR}/sorted_bravo_dbsnp_maf0.001_pass_vrt1.vcf.gz