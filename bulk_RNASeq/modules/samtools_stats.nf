process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'samtools_stats'

    input:
    tuple val(meta), path(input), path(input_index)
    path fasta

    output:
    tuple val(meta), path("*.stats"), emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input} \\
        > ${prefix}.stats
    """
}