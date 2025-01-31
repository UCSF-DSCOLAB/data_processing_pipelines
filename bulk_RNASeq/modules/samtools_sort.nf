process SAMTOOLS_SORT {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'samtools_sort', 'per_sample'
    memory {
        // File size in GB
        fileSize = bam.size() / (1024 * 1024 * 1024)
        return 3.GB + (1.GB * fileSize * 0.001)
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.csi"), emit: csi, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools sort $args -@ $task.cpus -o ${prefix}.bam -T $prefix $bam
    """
}
