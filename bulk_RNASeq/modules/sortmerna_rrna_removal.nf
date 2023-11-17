process SORTMERNA_RIBOSOMAL_RNA_REMOVAL {
    tag "$meta.id"
    label 'sortmerna_ribosomal_rna_removal'
    memory {
        if (meta.single_end) {
          // File size in GB
          fileSize = reads.size() / (1024 * 1024 * 1024)
        } else {
          // File size in GB
          fileSize = reads[0].size() / (1024 * 1024 * 1024)
        }
        if (fileSize > 3) {
            fileSize = 3
        }
        return 64.GB * (1 + (fileSize * 4))
    }

    input:
    tuple val(meta), path(reads)
    path  rrna_ref_fastas

    output:
    tuple val(meta), path("*.sortmerna.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def refs = "${rrna_ref_fastas.join(' --ref ')}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        sortmerna \\
            --ref $refs \\
            --reads $reads \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            $args

        mv non_rRNA_reads.f*q.gz ${prefix}.sortmerna.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log
        """
    } else {
        """
        sortmerna \\
            --ref $refs \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2 \\
            $args

        mv non_rRNA_reads_fwd.f*q.gz ${prefix}_R1.sortmerna.fastq.gz
        mv non_rRNA_reads_rev.f*q.gz ${prefix}_R2.sortmerna.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log
        """
    }
}