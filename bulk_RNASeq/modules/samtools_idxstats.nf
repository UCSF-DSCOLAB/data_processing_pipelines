process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'samtools_idxstats'

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