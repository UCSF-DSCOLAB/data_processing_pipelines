process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'samtools_idxstats', 'per_sample'
    memory {
        // File size in GB
        fileSize = bam.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.0001)
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        idxstats \\
        $bam \\
        > ${prefix}.idxstats
    """
}