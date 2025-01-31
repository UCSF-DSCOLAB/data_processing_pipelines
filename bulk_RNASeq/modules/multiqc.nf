process MULTIQC {
    tag "MultiQC"
    publishDir "${params.results_directory}/multiqc", mode: 'copy'
    // clusterOptions = '-S /bin/bash'
    label 'multiqc'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data

    script:
    def files = multiqc_files.collect { it.toString() }
    if (files.size() > 0) {
        """
        multiqc -f ${files.join(' ')} -o .
        """
    } else {
        echo "No files to process with MultiQC."
        """
        touch empty_multiqc_report.txt
        """
    }
}