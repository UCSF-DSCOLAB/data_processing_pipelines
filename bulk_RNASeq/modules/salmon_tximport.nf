process SALMON_TXIMPORT {
    publishDir "${params.results_directory}/salmon", mode: 'copy'
    clusterOptions = '-S /bin/bash'

    input:
    path ("salmon/*")
    path  tx2gene

    output:
    path "*gene_tpm.tsv"                 , emit: tpm_gene
    path "*gene_counts.tsv"              , emit: counts_gene
    path "*gene_counts_length_scaled.tsv", emit: counts_gene_length_scaled
    path "*gene_counts_scaled.tsv"       , emit: counts_gene_scaled
    path "*transcript_tpm.tsv"           , emit: tpm_transcript
    path "*transcript_counts.tsv"        , emit: counts_transcript

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tximport.r \\
        NULL \\
        salmon \\
        salmon.merged
    """
}