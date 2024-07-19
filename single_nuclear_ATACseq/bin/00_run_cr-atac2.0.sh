#!/bin/bash
set -e
set -o nounset

mkdir ${TMPDIR}/cellranger_atac_count_${SAMPLE}
cd ${TMPDIR}/cellranger_atac_count_${SAMPLE}

if [[ $(grep -c "ATAC Seq" ${LIBRARIES_CSV}) -ne 1 ]]
then
  echo "ERROR: ${LIBRARIES_CSV} must contain ONLY one library of type 'ATAC Seq'"
  exit 1
fi

fastq_dir=$(grep "ATAC Seq" ${LIBRARIES_CSV} | cut -f1 -d ",")

singularity exec \
            -B ${fastq_dir} \
            -B ${REFERENCE} \
            -B ${PWD} \
            --pwd ${PWD} \
            /krummellab/data1/singularity_images/cellranger-atac/${CELLRANGERATACVERSION}/cellranger-atac.sif \
                cellranger-atac-${CELLRANGERATACVERSION} count \
                                                --id=${SAMPLE} \
                                                --fastqs=${fastq_dir} \
                                                --sample=${SAMPLE} \
                                                --reference=${REFERENCE} \
                                                --localcores=32 \
                                                --localmem=150

if [ -f ${OUTDIR} ]
then
    mv ${TMPDIR}/cellranger_atac_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_counts
else
    mkdir -p $(dirname ${OUTDIR})
    #mv ${TMPDIR}/cellranger_atac_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
    mv ${TMPDIR}/cellranger_atac_count_${SAMPLE}/  ${OUTDIR}
fi

