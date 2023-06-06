process KALLISTO_QUANT {
    tag "$meta.id"
    publishDir "${params.results_directory}/kallisto", mode: 'copy'
    cpus 12
    memory '64 GB'
    conda "$baseDir/envs/kallisto.yml"

    input:
    tuple val(meta), path(reads)

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
        kallisto quant --index ${params.genome_dir}/kallisto --output-dir ${prefix} -t ${task.cpus} --single -l ${params.fragment_length_mean} -s ${params.fragment_length_std} ${reads} > ${prefix}.kallisto.log
        mv ${prefix}/abundance.h5 ${prefix}.kallisto.abundance.h5
        mv ${prefix}/abundance.tsv ${prefix}.kallisto.abundance.tsv
        mv ${prefix}/run_info.json ${prefix}.kallisto.run_info.json
        """
    } else {
        """
        kallisto quant --index ${params.genome_dir}/kallisto --output-dir ${prefix} -t ${task.cpus} ${reads[0]} ${reads[1]} > ${prefix}.kallisto.log
        mv ${prefix}/abundance.h5 ${prefix}.kallisto.abundance.h5
        mv ${prefix}/abundance.tsv ${prefix}.kallisto.abundance.tsv
        mv ${prefix}/run_info.json ${prefix}.kallisto.run_info.json
        """
    }
}