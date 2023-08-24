process CUSTOM_MERGE_COUNTS {
    tag "$samplesheet"
    label 'custom_merge_counts'
    publishDir "${params.results_directory}/merged_results", mode: 'copy'

    input:
    path counts

    output:
    path 'merged_counts.tsv'       , emit: merged_counts

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnavar/bin/
    """
    merge_kallisto_counts.py \\
        ${counts.join(',')} \\
        merged_counts.tsv
    """
}
