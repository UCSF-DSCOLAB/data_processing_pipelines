process EXTRACT_MAPPED_READS {
    tag "$meta.id"
    cpus 32
    memory '64 GB'
    conda "$baseDir/envs/samtools.yml"

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