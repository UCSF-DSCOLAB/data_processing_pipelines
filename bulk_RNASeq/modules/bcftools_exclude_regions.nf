process BCFTOOLS_EXCLUDE_REGIONS {
    tag "cohort"
    label 'bcftools_demux_filter'
    publishDir "${params.results_directory}/snps/demux_ready", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path excluded_regions

    output:
    tuple val(meta), path("cohort.demux_ready.final.vcf.gz"), path("cohort.demux_ready.final.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view \\
        --targets-file ^$excluded_regions \\
        --output-type z \\
        --output cohort.demux_ready.final.vcf.gz \\
        $vcf
    bcftools index --tbi cohort.demux_ready.final.vcf.gz
    """
}
