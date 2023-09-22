process SALMON_TX2GENE {
    tag "$gtf"
    publishDir "${params.results_directory}/salmon", mode: 'copy'

    input:
    path ("salmon/*")
    path gtf

    output:
    path "*.tsv"       , emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tx2gene.py \\
        --gtf $gtf \\
        --salmon salmon \\
        --id $params.gtf_group_features \\
        --extra $params.gtf_extra_attributes \\
        -o salmon_tx2gene.tsv
    """
}