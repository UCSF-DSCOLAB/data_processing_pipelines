process GATK4_BASE_RECALIBRATOR {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'gatk4_recalibrator', 'per_sample'
    memory {
        // File size in GB
        fileSize = input.size() / (1024 * 1024 * 1024)
        return 5.GB + (1.GB * fileSize * 0.1)
    }

    input:
    tuple val(meta), path(input)
    tuple val(meta), path(bai)
    path  fasta
    path  fai
    path  dict
    path  known_sites
    path  known_sites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def interval_command = intervals ? "--intervals $intervals" : ""
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" BaseRecalibrator  \\
        --input $input \\
        --output ${prefix}.table \\
        --reference $fasta \\
        $sites_command \\
        --tmp-dir \$PWD \\
        $args
    """
}