process GATK4_SPLITNCIGARREADS {
    tag "$meta.id"
    label 'gatk4_splitncigarreads'

    input:
    tuple val(meta), path(bam), path(bai)
    path  genome
    path  genome_idx
    path  dict

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
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SplitNCigarReads \\
        --input $bam \\
        --output ${prefix}.bam \\
        --reference $genome \\
        --tmp-dir $params.tmp_dir \\
        $args
    """
}