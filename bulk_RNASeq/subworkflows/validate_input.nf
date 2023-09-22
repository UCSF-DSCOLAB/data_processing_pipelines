nextflow.enable.dsl=2
//
// Check input samplesheet and get read channels
//

include { VALIDATE_SAMPLE_SHEET } from '../modules/validate_sample_sheet'

// params.sample_sheet = '/krummellab/data1/alaa/data/tests/bulk_rnaseq/pipe_dev/test_sample_sheet.csv'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    VALIDATE_SAMPLE_SHEET (
        samplesheet
    )
    .csv
    .splitCsv ( header:true, sep:',' )
    .map { create_fastq_channel(it) }
    .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
