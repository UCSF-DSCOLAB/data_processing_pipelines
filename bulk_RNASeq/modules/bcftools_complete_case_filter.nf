process BCFTOOLS_COMPLETE_CASE_FILTER {
    tag "cohort"
    label 'bcftools_demux_filter'
    publishDir "${params.results_directory}/snps/demux_ready", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("cohort.demux_ready.complete.vcf.gz"), path("cohort.demux_ready.complete.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools +fill-tags $vcf -Ou -- -t F_MISSING | \\
        bcftools view \\
            --include 'F_MISSING=0' \\
            --output-type z \\
            --output cohort.demux_ready.complete.vcf.gz
    bcftools index --tbi cohort.demux_ready.complete.vcf.gz
    """
}
