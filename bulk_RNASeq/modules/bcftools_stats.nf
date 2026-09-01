process BCFTOOLS_STATS {
    tag "${label}"
    label 'bcftools_stats'
    publishDir "${params.results_directory}/snps/qc", mode: 'copy'

    input:
    tuple val(label), path(vcf), path(tbi)

    output:
    tuple val(label), path("*.bcftools.stats.txt"), emit: stats

    script:
    """
    bcftools stats --samples - $vcf > ${label}.bcftools.stats.txt
    """
}
