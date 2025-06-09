process CUSTOM_MERGE_COUNTS {
    tag "merge_counts"
    // clusterOptions = '-S /bin/bash'
    label 'custom_merge_counts'
    publishDir "${params.results_directory}/merged_results", mode: 'copy'

    input:
    path counts_files

    output:
    path 'merged_counts.tsv', emit: merged_counts

    script:
    """
    merge_kallisto_counts.py \\
        ${counts_files.join(' ')} \\
        merged_counts.tsv
    """
}