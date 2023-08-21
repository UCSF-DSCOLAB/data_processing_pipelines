process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'picard_markduplicates'
    cpus 2
    memory '31 GB'
    conda "$baseDir/envs/gatk.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , optional:true, emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk \\
        --java-options "-Xmx${task.memory.toGiga()}g" MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.picard.bam \\
        --REFERENCE_SEQUENCE $fasta \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt
    """
}
