#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=1:00:00
#SBATCH -e /krummellab/data1/%u/logs/run_freemux-%j.err
#SBATCH -o /krummellab/data1/%u/logs/run_freemux-%j.out


LIBRARY=$1
NUMSAMPLES=$2

CONTAINER="/krummellab/data1/singularity_images/popscle/da70fc78da385ef049e0e890342acfd62842cae0/popscle.sif"

SAMPLE=${LIBRARY}
BARCODE_DIR=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/cell_filter/
BARCODES=${BARCODE_DIR}/post_amulet_barcodes_of_interest_filt.list

FMX_IN=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet/
FMX_OUT=/krummellab/data1/immunox/${PROJECT}/data/single_nuclear_ATAC/processed/${LIBRARY}/freemuxlet_n${NUMSAMPLES}/
mkdir -p ${FMX_OUT}
singularity exec \
                -B ${FMX_IN} -B ${FMX_OUT} -B ${BARCODE_DIR} \
                ${CONTAINER} popscle freemuxlet \
                        --plp ${FMX_IN}/${SAMPLE} \
                        --out ${FMX_OUT}/${SAMPLE} \
                        --nsample ${NUMSAMPLES} \
                        --seed 212 \
                        --group-list ${BARCODES}
