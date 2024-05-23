#!/bin/bash

POOL=$1
REF_VCF_FILE=$2
FMX_VCF_FILE=$3
REF_TYPE=$4

cp ${REF_VCF_FILE} ground_truth.vcf.gz
bcftools view -O b -o ground_truth.bcf ground_truth.vcf.gz
bcftools view -O b -o ${POOL}.bcf ${FMX_VCF_FILE}
bcftools index ground_truth.bcf
bcftools index ${POOL}.bcf

bcftools filter --exclude "GQ<20" -O b -o ${POOL}_filtered.bcf ${POOL}.bcf
bcftools index ${POOL}_filtered.bcf

if [[ "${REF_TYPE}" == "array" ]]
then
    bcftools gtcheck -g ground_truth.bcf ${POOL}_filtered.bcf | grep -v "^INFO" > ${POOL}_gtcheck.out
else
    bcftools isec -O b -p ${POOL}_with_ground_truth -n =2 ${POOL}_filtered.bcf ground_truth.bcf
    bcftools merge --merge snps -o pool_with_ground_truth.bcf -O b ${POOL}_with_ground_truth/0000.bcf \
           ${POOL}_with_ground_truth/0001.bcf
    bcftools index pool_with_ground_truth.bcf
    bcftools gtcheck pool_with_ground_truth.bcf | grep -v "^INFO" > ${POOL}_gtcheck.out
fi