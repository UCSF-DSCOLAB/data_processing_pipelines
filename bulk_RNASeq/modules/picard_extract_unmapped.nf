process EXTRACT_UNMAPPED_READS {
    tag "$meta.id"
    cpus 32
    memory '64 GB'
    publishDir "${params.results_directory}/star", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.trimmed.unmapped.{1,2}.fastq.gz'), emit: unmapped_fastq


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}."
    """
    java -Xmx${task.memory.toGiga()-5}g \
        -jar /opt/picard/picard.jar \
        ViewSam \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNMENT_STATUS=Unaligned \
        PF_STATUS=All \
        I=${bam} \
        > ${prefix}.trimmed.unmapped.bam

    java -Xmx${task.memory.toGiga()-5}g \
        -jar /opt/picard/picard.jar \
        SamToFastq \
        VALIDATION_STRINGENCY=SILENT \
        I=${prefix}.trimmed.unmapped.bam \
        FASTQ=${prefix}.trimmed.unmapped.1.fq.gz \
        SECOND_END_FASTQ=${prefix}.trimmed.unmapped.2.fq.gz
    """
}