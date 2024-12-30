process BCFTOOLS_CONTIG_CONVERSION {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    label 'bcftools_contig_conversion', 'per_sample'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcf.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize)
    }

    input:
    tuple val(meta), path(vcf)
    path format_map

    output:
    tuple val(meta), path("*.formatted.vcf.gz"), emit: formatted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools annotate \\
            --rename-chrs $format_map \\
            --threads $task.cpus -Oz \\
            --output ${prefix}.formatted.vcf.gz \\
            $vcf
            $args
    """
}