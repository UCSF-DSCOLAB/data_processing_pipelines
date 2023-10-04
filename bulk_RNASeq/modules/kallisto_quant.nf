process KALLISTO_QUANT {
    tag "$meta.id"
    label 'kallisto_quant'
    publishDir "${params.results_directory}/kallisto", mode: 'copy'
    memory {
        if {meta.single_end} {
            // File size in GB
            fileSize = reads.size() / (1024 * 1024 * 1024)
        } else {
            // File size in GB
            fileSize = reads[0].size() / (1024 * 1024 * 1024)
        }
        return 7.GB * (1 + (fileSize * 0.5))
    }

    input:
    tuple val(meta), path(reads)
    path transcript_index

    output:
    tuple val(meta), path("${prefix}.kallisto.abundance.h5") , emit: abundance_h5
    tuple val(meta), path("${prefix}.kallisto.abundance.tsv") , emit: abundance_tsv
    tuple val(meta), path("${prefix}.kallisto.run_info.json"), emit: json_info, optional: true
    tuple val(meta), path("${prefix}.kallisto.log"), emit: log, optional: true
    
    
    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        kallisto quant --index $transcript_index --output-dir ${prefix} -t ${task.cpus} --single -l ${params.fragment_length_mean} -s ${params.fragment_length_std} ${reads} > ${prefix}.kallisto.log
        mv ${prefix}/abundance.h5 ${prefix}.kallisto.abundance.h5
        mv ${prefix}/abundance.tsv ${prefix}.kallisto.abundance.tsv
        mv ${prefix}/run_info.json ${prefix}.kallisto.run_info.json
        """
    } else {
        """
        kallisto quant --index $transcript_index --output-dir ${prefix} -t ${task.cpus} ${reads[0]} ${reads[1]} > ${prefix}.kallisto.log
        mv ${prefix}/abundance.h5 ${prefix}.kallisto.abundance.h5
        mv ${prefix}/abundance.tsv ${prefix}.kallisto.abundance.tsv
        mv ${prefix}/run_info.json ${prefix}.kallisto.run_info.json
        """
    }
}