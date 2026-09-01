process BCFTOOLS_DEMUX_FILTER {
    tag "cohort"
    label 'bcftools_demux_filter'
    publishDir "${params.results_directory}/snps/demux_ready", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai

    output:
    tuple val(meta), path("cohort.demux_ready.vcf.gz"), path("cohort.demux_ready.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools norm \\
        --fasta-ref $fasta \\
        --multiallelics -any \\
        --output-type z \\
        --output cohort.normalized.vcf.gz \\
        $vcf
    bcftools index --tbi cohort.normalized.vcf.gz
    bcftools +fill-tags cohort.normalized.vcf.gz -Ou -- -t AC,AN,AF,F_MISSING | \\
        bcftools view \\
            --include 'F_MISSING <= ${params.demux_max_missing_fraction} && AC > 0 && AC < AN' \\
            --output-type z \\
            --output cohort.demux_ready.vcf.gz
    bcftools index --tbi cohort.demux_ready.vcf.gz
    """
}
