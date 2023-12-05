process GATK4_APPLY_BQSR {
    tag "$meta.id"
    label 'gatk4_apply_bqsr'
    memory {
        // File size in GB
        fileSize = input.size() / (1024 * 1024 * 1024)
        return 1.GB + (2.GB * fileSize * 0.1)
    }

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table)
    path  genome
    path  genome_idx
    path  genome_dict

    output:
    tuple val(meta), path("*_bqsr.bam") , emit: bam
    // tuple val(meta), path("*.bai") , emit: bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyBQSR \\
        --input $input \\
        --output ${prefix}_bqsr.bam \\
        --reference $genome \\
        --bqsr-recal-file $bqsr_table \\
        --tmp-dir \$PWD \\
        $args
    """
}