process VALIDATE_SAMPLE_SHEET {
    tag "$samplesheet"

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
