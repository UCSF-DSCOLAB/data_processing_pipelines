process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'picard_markduplicates', 'per_sample'
    memory {
        // File size in GB
        fileSize = bam.size() / (1024 * 1024 * 1024)
        return 20.GB + (5.GB * fileSize * 4)
    }

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
