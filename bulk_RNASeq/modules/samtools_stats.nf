process SAMTOOLS_STATS {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    label 'samtools_stats', 'per_sample'
    memory {
        // File size in GB
        fileSize = input.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.0001)
    }

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