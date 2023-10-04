process FASTP_TRIM_ADAPTERS {
    tag "$meta.id"
    label 'fastp_trim_adapters'
    memory {
        // File size in GB
        fileSize = reads.size() / (1024 * 1024 * 1024)
        return 1.GB + (1.GB * fileSize)
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fastq.gz") , emit: trimmed_reads
    path("*.fastp.json")                      , emit: json_report
    path("*.fastp.html")                      , emit: html_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_sequence_1 = params.adapter_sequence_1 ? "--adapter_sequence ${params.adapter_sequence_1}" : ""
    def adapter_sequence_2 = params.adapter_sequence_2 ? "--adapter_sequence_r2 ${params.adapter_sequence_2}" : ""
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
        fastp \\
          --in1 ${prefix}.fastq.gz \\
          --out1 ${prefix}.fastp.fastq.gz \\
          --length_required 20 \
          $adapter_sequence_1 \\
          $adapter_sequence_2 \\
          --correction  \
          --trim_poly_g  \
          --thread ${task.cpus} \\
          --json ${prefix}.fastp.json \\
          --html ${prefix}.fastp.html
        """
    } else {
        """
        [ ! -f  ${prefix}_R1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f  ${prefix}_R2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_R2.fastq.gz
        fastp \\
          --in1 ${prefix}_R1.fastq.gz \\
          --in2 ${prefix}_R2.fastq.gz \\
          --out1 ${prefix}_R1.fastp.fastq.gz \\
          --out2 ${prefix}_R2.fastp.fastq.gz \\
          --length_required 20 \
          $adapter_sequence_1 \\
          $adapter_sequence_2 \\
          --correction  \
          --trim_poly_g  \
          --thread ${task.cpus} \\
          --json ${prefix}.fastp.json \\
          --html ${prefix}.fastp.html
        """
    }
}