process GATK4_SPLITNCIGARREADS {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'gatk4_splitncigarreads', 'per_sample'
    memory {
        // File size in GB
        fileSize = bam.size() / (1024 * 1024 * 1024)
        return 200.GB + (1.GB * fileSize * 5)
    }

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
        --tmp-dir \$PWD \\
        $args
    """
}