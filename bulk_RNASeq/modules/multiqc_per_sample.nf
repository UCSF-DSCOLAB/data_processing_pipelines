process MULTIQC_PER_SAMPLE {
    tag "$meta.id"
    publishDir "${params.results_directory}/multiqc_per_sample", mode: 'copy'
    label 'multiqc_per_sample'
    memory {
        // File size in GB
        fileSize = log_files.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize)
    }

    input:
    tuple val(meta), path(log_files)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    multiqc --filename ${prefix}_multiqc_report.html -f $args .
    """
}