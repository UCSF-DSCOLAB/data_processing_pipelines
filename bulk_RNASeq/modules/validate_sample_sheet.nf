process VALIDATE_SAMPLE_SHEET {
    tag "${samplsheet}"
    label 'validate_sample_sheet'
    clusterOptions = "-S /bin/bash"
    memory {
        // File size in GB
        fileSize = samplesheet.size() / (1024 * 1024 * 1024)
        return 5.GB * (1+fileSize)
    }

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnavar/bin/
    """
    validate_sample_sheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """
}
