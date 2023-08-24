process SAMTOOLS_BAM_TO_CRAM {
    tag "$meta.id"
    label 'samtools_bam_to_cram'
    publishDir "${params.results_directory}/star", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path reference_fasta
    val prefix_addon

    output:
    tuple val(meta), path('*.mapped.cram'), emit: mapped_cram

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -@ ${task.cpus} \\
        -C \\
        --no-PG \\
        -T ${reference_fasta} \\
        -o ${prefix}${prefix_addon}.mapped.cram \\
        ${bam}
    """
}