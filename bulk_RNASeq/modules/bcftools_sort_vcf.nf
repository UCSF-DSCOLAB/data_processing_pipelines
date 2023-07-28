process BCFTOOLS_SORT_VCF {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    conda "$baseDir/envs/bcftools.yml"

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
    bcftools sort $vcf -o ${prefix}.sorted.vcf.gz -O z
    """
}