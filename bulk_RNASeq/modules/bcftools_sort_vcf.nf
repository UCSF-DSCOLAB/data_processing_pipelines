process BCFTOOLS_SORT_VCF {
    tag "$meta.id"
    label 'bcftools_sort_vcf'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcf.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.01)
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.sorted.vcf.gz"), emit: sorted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools sort \\
            --output ${prefix}.sorted.vcf.gz -Oz \\
            $vcf
    """
}