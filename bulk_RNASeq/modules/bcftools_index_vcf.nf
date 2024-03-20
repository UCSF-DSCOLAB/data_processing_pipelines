process BCFTOOLS_INDEX_VCF {
    tag "$meta.id"
    label 'bcftools_index_vcf'
    scratch = false
    publishDir "${params.results_directory}/snps", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcf.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize * 0.01)
    }

    input:
    tuple val(meta), path(vcf)

    output:
    // tuple val(meta), path("*.formatted.vcf.gz"), emit: sorted_vcf
    tuple val(meta), path("*.sorted.vcf.gz.tbi"), emit: vcf_index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools index --tbi $vcf
    """
}