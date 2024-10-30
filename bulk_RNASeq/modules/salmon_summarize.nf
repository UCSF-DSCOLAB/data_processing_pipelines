process SALMON_SUMMARIZE_EXPERIMENT {
    tag "$tx2gene"
    clusterOptions = '-S /bin/bash'
    publishDir "${params.results_directory}/salmon", mode: 'copy'

    input:
    path counts
    path tpm
    path tx2gene

    output:
    path "*.rds"       , emit: rds

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_summarize_experiment.r \\
        NULL \\
        $counts \\
        $tpm
    """
}
