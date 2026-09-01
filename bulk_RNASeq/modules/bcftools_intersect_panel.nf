process BCFTOOLS_INTERSECT_PANEL {
    tag "cohort"
    label 'bcftools_demux_filter'
    publishDir "${params.results_directory}/snps/demux_ready", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path panel
    path panel_tbi
    path fasta
    path fai

    output:
    tuple val(meta), path("cohort.demux_ready.panel.vcf.gz"), path("cohort.demux_ready.panel.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools norm \\
        --check-ref e \\
        --fasta-ref $fasta \\
        --multiallelics -any \\
        --output-type z \\
        --output normalized.panel.vcf.gz \\
        $panel
    bcftools index --tbi normalized.panel.vcf.gz
    bcftools isec \\
        --collapse none \\
        --nfiles 2 \\
        --write 1 \\
        --output-type z \\
        --output cohort.demux_ready.panel.vcf.gz \\
        $vcf normalized.panel.vcf.gz
    bcftools index --tbi cohort.demux_ready.panel.vcf.gz
    """
}
