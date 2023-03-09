process GATK4_SPLITNCIGARREADS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path  fasta
    path  fai
    path  dict
    path  tmpDir

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.bai'), emit: bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def interval_command = intervals ? "--intervals $intervals" : ""
    """
    gatk --java-options "-Xmx314g" SplitNCigarReads \\
        --input $bam \\
        --output ${prefix}.bam \\
        --reference $fasta \\
        $args
    """
}