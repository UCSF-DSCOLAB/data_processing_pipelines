process GENOTYPE_VCF_QC {
    tag "cohort"
    label 'genotype_vcf_qc'
    publishDir "${params.results_directory}/snps/qc", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    path "cohort.sample_genotype_qc.tsv", emit: sample_qc
    path "cohort.pairwise_discrimination.tsv", emit: pairwise_qc

    script:
    """
    summarize_genotype_vcf.py \\
        $vcf \\
        --sample-output cohort.sample_genotype_qc.tsv \\
        --pairwise-output cohort.pairwise_discrimination.tsv \\
        --minimum-discordant ${params.demux_min_pairwise_discordant_snps}
    """
}
