#!/bin/bash
#
# Copyright (C): 2020 - Gert Hulselmans
#
# Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping only reads:
#            - which overlap with SNPs in the VCF file
#            - and which have a cell barcode contained in the cell barcode list
#          Keeping only relevant reads for dsc-pileup can speedup it up several hunderd times.



# Function to check if any of the programs in a pipe failed.
check_exit_codes () {
    local GET_PIPESTATUS="${PIPESTATUS[@]}";
    local exit_code;

    for exit_code in ${GET_PIPESTATUS} ; do
        if [ ${exit_code} -ne 0 ] ; then
             return ${exit_code};
        fi
    done

    return 0;
}


filter_bam_file_for_popscle_dsc_pileup () {
    local input_bam_filename="${1}";
    local barcodes_tsv_filename="${2}";
    local vcf_filename="${3}";
    local output_bam_filename="${4}";

    local exit_code=0;

    if [ ${#@} -ne 4 ] ; then
        printf 'Usage:   filter_bam_file_for_popscle_dsc_pileup input_bam_filename barcodes_tsv_filename vcf_filename output_bam_filename\n\n';
        printf 'Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping reads:\n';
        printf '           - which overlap with SNPs in the VCF file\n';
        printf '           - and which have a cell barcode contained in the cell barcode list\n';
        printf '         Keeping only relevant reads for popscle dsc-pileup can speedup it up several hunderd times.\n\n';

        return 1;
    fi

    if [ ! -f  "${input_bam_filename}" ] ; then
        printf 'Error: Input (CellRanger) BAM file "%s" could not be found.\n' "${input_bam_filename}";
        return 2;
    fi

    if [ ! -f  "${barcodes_tsv_filename}" ] ; then
        printf 'Error: File with barcodes "%s" could not be found.\n' "${barcodes_tsv_filename}";
        return 2;
    fi

    if [ ! -f  "${vcf_filename}" ] ; then
        printf 'Error: File with unique SNPs per sample "%s" could not be found.\n' "${vcf_filename}";
        return 2;
    fi


    # Create much smaller BAM file for dsc-pileup of popscle:
    #   - Convert VCF file with unique SNPs for each sample
    #     to a BED file and merge adjacent SNP regions to one.
    #   - Only include reads that contain a SNP position
    #     and which contain a cell barcode of interest.
    if [ "${barcodes_tsv_filename%.gz}".gz = "${barcodes_tsv_filename}" ] ; then
        # Barcodes file is compressed with gzip.
        bedtools merge -i "${vcf_filename}" \
          | samtools view\
                -@ 8 \
                -L - \
                -D CB:<(zcat "${barcodes_tsv_filename}") \
                -o "${output_bam_filename}" \
                "${input_bam_filename}";

        # Check if any of the previous commands failed.
        check_exit_codes;

        exit_code=$?;
    else
        # Barcodes file is uncompressed.
        bedtools merge -i "${vcf_filename}" \
          | samtools view\
                -@ 8 \
                -L - \
                -D CB:"${barcodes_tsv_filename}" \
                -o "${output_bam_filename}" \
                "${input_bam_filename}";

        # Check if any of the previous commands failed.
        check_exit_codes;

        exit_code=$?;
    fi


    return ${exit_code};
}

filter_bam_file_for_popscle_dsc_pileup "${@}";