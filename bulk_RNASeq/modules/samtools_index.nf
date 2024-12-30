process SAMTOOLS_INDEX {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    label 'samtools_index', 'per_sample'
    memory {
        // File size in GB
        fileSize = input.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.0001)
    }

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input
    """
}