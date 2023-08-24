process PICARD_EXTRACT_UNMAPPED_READS {
    tag "$meta.id"
    label 'picard_extract_unmapped_reads'
    publishDir "${params.results_directory}/star", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.trimmed.unmapped.{1,2}.fastq.gz'), emit: unmapped_fastq


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    picard ViewSam \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNMENT_STATUS=Unaligned \
        PF_STATUS=All \
        I=${bam} \
        > ${prefix}.trimmed.unmapped.bam

    picard SamToFastq \
        VALIDATION_STRINGENCY=SILENT \
        I=${prefix}.trimmed.unmapped.bam \
        FASTQ=${prefix}.trimmed.unmapped.1.fq.gz \
        SECOND_END_FASTQ=${prefix}.trimmed.unmapped.2.fq.gz
    """
}