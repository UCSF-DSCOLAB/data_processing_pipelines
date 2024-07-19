#!/bin/bash
REF_VCF=/krummellab/data1/immunox/${PROJECT}/genotypes/glands_genotyping/merged_all.vcf.gz
for LIB_PATH in `ls -d /krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/*`;
do 
  LIBRARY="${LIB_PATH##*/}"
  echo $LIBRARY
  bcftools view -O b -o ${LIBRARY}.clust1.bcf /krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet/${LIBRARY}.clust1.vcf.gz
  bcftools index ${LIBRARY}.clust1.bcf
  bcftools filter --exclude "GQ<20" -O b -o ${LIBRARY}.clust1_filtered.bcf ${LIBRARY}.clust1.bcf
  bcftools index ${LIBRARY}.clust1_filtered.bcf
  
  bcftools isec -O b -p ${LIBRARY}_with_ground_truth -n =2 ${LIBRARY}.clust1_filtered.bcf ${REF_VCF}
  
  cd ${LIBRARY}_with_ground_truth
  bcftools merge --merge snps -o pool_with_ground_truth.bcf -O b 0000.bcf 0001.bcf
  bcftools index pool_with_ground_truth.bcf
  
  bcftools gtcheck -G 1 pool_with_ground_truth.bcf > ../${LIBRARY}_gtcheck.out
  cd ..
  rm -rf ${LIBRARY}_with_ground_truth 
  mv ${LIBRARY}_gtcheck.out /krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet/
done