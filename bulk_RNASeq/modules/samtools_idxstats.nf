process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    conda "$baseDir/envs/samtools.yml"

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