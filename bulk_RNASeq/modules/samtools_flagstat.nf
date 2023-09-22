process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    label 'samtools_flagstat'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        $bam \\
        > ${prefix}.flagstat
    """
}
