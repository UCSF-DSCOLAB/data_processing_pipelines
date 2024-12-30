// Merge fastq files, if needed
process CAT_FASTQ {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    label 'cat_fastq', 'per_sample'

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads.collect{ it.toString() }
    if (meta.single_end) {
        if (readList.size > 1) {
            """
            cat ${readList.join(' ')} > ${prefix}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${prefix}_2.merged.fastq.gz
            """
        }
    }
}