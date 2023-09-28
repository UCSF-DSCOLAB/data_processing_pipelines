#!/bin/bash

POOL=$1
REF_VCF_FILE=$2
FMX_VCF_FILE=$3

REF_VCF=${REF_VCF_FILE%%.*}
FMX_VCF=${FMX_VCF_FILE%%.*}

cp ${REF_VCF_FILE} ground_truth.vcf.gz
bcftools view -O b -o ground_truth.bcf ground_truth.vcf.gz
bcftools view -O b -o ${FMX_VCF}.bcf ${FMX_VCF_FILE}
bcftools index ground_truth.bcf
bcftools index ${FMX_VCF}.bcf

bcftools filter --exclude "GQ<20" -O b -o ${FMX_VCF}_filtered.bcf ${FMX_VCF}.bcf
bcftools index ${FMX_VCF}_filtered.bcf

bcftools isec -O b -p ${POOL}_with_ground_truth -n =2 ${FMX_VCF}_filtered.bcf ground_truth.bcf

bcftools merge --merge snps -o pool_with_ground_truth.bcf -O b ${POOL}_with_ground_truth/0000.bcf \
    ${POOL}_with_ground_truth/0001.bcf
bcftools index pool_with_ground_truth.bcf

bcftools gtcheck -G 1 pool_with_ground_truth.bcf > ${POOL}_gtcheck.out
