process GATK4_APPLY_BQSR {
    tag "$meta.id"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table)
    path  genome
    path  genome_idx
    path  genome_dict

    output:
    tuple val(meta), path("*.bam") , emit: bam
    tuple val(meta), path("*.bai") , emit: bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk --java-options "-Xmx314g" ApplyBQSR \\
        --input $input \\
        --output ${prefix}.bam \\
        --reference $genome \\
        --bqsr-recal-file $bqsr_table \\
        --tmp-dir . \\
        $args
    # Index the BAM file
    samtools index ${prefix}.bam
    """
}