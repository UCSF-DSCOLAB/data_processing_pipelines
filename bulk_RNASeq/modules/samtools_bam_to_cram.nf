process SAMTOOLS_BAM_TO_CRAM {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'samtools_bam_to_cram'
    publishDir "${params.results_directory}/star", mode: 'copy'
    memory {
        def sizeGiB   = Math.ceil( bam.size() / (1024 ** 3) )
        def required  = 1.GB + sizeGiB * 1.GB
        return required < 8.GB ? 8.GB : required   // equivalent to Math.max(required, 8.GB)
    }

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