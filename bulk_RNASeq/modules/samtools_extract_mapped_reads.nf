process SAMTOOLS_EXTRACT_MAPPED_READS {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'samtools_extract_mapped_reads', 'per_sample'
    publishDir "${params.results_directory}/star", mode: 'copy'
    memory {
        // File size in GB
        fileSize = bam.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.001)
    }

    input:
    tuple val(meta), path(bam)
    val prefix_addon

    output:
    tuple val(meta), path('*.mapped.bam'), emit: mapped_bam

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view  -b \\
        -F 0x4 \\
        -@ ${task.cpus} \\
        --no-PG \\
        -o ${prefix}${prefix_addon}.mapped.bam \\
        ${bam}
    """
}