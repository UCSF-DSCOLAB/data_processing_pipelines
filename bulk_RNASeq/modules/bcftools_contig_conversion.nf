process BCFTOOLS_CONTIG_CONVERSION {
    tag "$meta.id"
    label 'bcftools_contig_conversion'
    scratch = false

    publishDir "${params.results_directory}/snps", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcf.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.01)
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