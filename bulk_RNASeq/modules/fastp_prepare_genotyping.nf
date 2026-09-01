process FASTP_PREPARE_GENOTYPING {
    tag "$meta.id"
    label 'fastp_trim_adapters', 'per_sample'
    publishDir "${params.results_directory}/genotyping_reads", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.genotyping.fastp.fastq.gz"), emit: trimmed_reads
    path("*.genotyping.fastp.json"), emit: json_report
    path("*.genotyping.fastp.html"), emit: html_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_sequence_1 = params.adapter_sequence_1 ? "--adapter_sequence ${params.adapter_sequence_1}" : ""
    def adapter_sequence_2 = params.adapter_sequence_2 ? "--adapter_sequence_r2 ${params.adapter_sequence_2}" : ""
    if (meta.single_end) {
        """
        fastp \\
          --in1 $reads \\
          --out1 ${prefix}.genotyping.fastp.fastq.gz \\
          --length_required 20 \\
          $adapter_sequence_1 \\
          --trim_poly_g \\
          --thread ${task.cpus} \\
          --json ${prefix}.genotyping.fastp.json \\
          --html ${prefix}.genotyping.fastp.html
        """
    } else {
        """
        fastp \\
          --in1 ${reads[0]} \\
          --in2 ${reads[1]} \\
          --out1 ${prefix}_R1.genotyping.fastp.fastq.gz \\
          --out2 ${prefix}_R2.genotyping.fastp.fastq.gz \\
          --length_required 20 \\
          $adapter_sequence_1 \\
          $adapter_sequence_2 \\
          --trim_poly_g \\
          --thread ${task.cpus} \\
          --json ${prefix}.genotyping.fastp.json \\
          --html ${prefix}.genotyping.fastp.html
        """
    }
}
