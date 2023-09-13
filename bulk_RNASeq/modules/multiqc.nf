process MULTIQC {
    publishDir "${params.results_directory}/multiqc", mode: 'copy'
    label 'multiqc'    

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f $args .
    """
}