process SORTMERNA_RIBOSOMAL_RNA_REMOVAL {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'sortmerna_ribosomal_rna_removal', 'per_sample'
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
        return 15.GB * (1 + (fileSize * 0.1))
    }
    publishDir "${params.results_directory}/trimmed_cleaned_reads", mode: 'copy'

    containerOptions "-B /scratch/"

    input:
    tuple val(meta), path(reads)
    path  rrna_ref_fastas

    output:
    tuple val(meta), path("*.sortmerna.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def refs   = "${rrna_ref_fastas.join(' --ref ')}"
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        set -euo pipefail

        count_reads () {
          local f="\$1"
          # Count FASTQ records (lines/4), works for .gz and plain
          if [[ "\$f" == *.gz ]]; then
            local l
            l=\$(zcat -f -- "\$f" 2>/dev/null | wc -l || echo 0)
            echo \$(( l / 4 ))
          else
            local l
            l=\$(wc -l < "\$f" 2>/dev/null || echo 0)
            echo \$(( l / 4 ))
          fi
        }

        R=\$(count_reads "$reads")
        if [[ "\$R" -eq 0 ]]; then
          echo "[sortmerna] No reads detected in $reads — skipping rRNA removal." > "${prefix}.sortmerna.log"
          # Emit an empty gzipped FASTQ to satisfy downstream expectations
          gzip -cn /dev/null > "${prefix}.sortmerna.fastq.gz"
          exit 0
        fi

        sortmerna \\
            --ref $refs \\
            --reads $reads \\
            --threads $task.cpus \\
            --workdir \$TMPDIR/ \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            $args

        mv non_rRNA_reads.f*q.gz "${prefix}.sortmerna.fastq.gz"
        mv rRNA_reads.log        "${prefix}.sortmerna.log"
        """
    } else {
        """
        set -euo pipefail

        count_reads () {
          local f="\$1"
          if [[ "\$f" == *.gz ]]; then
            local l
            l=\$(zcat -f -- "\$f" 2>/dev/null | wc -l || echo 0)
            echo \$(( l / 4 ))
          else
            local l
            l=\$(wc -l < "\$f" 2>/dev/null || echo 0)
            echo \$(( l / 4 ))
          fi
        }

        R1=\$(count_reads "${reads[0]}")
        R2=\$(count_reads "${reads[1]}")

        if [[ "\$R1" -eq 0 || "\$R2" -eq 0 ]]; then
          echo "[sortmerna] One or both mates have zero reads (R1=\$R1, R2=\$R2) — skipping rRNA removal." > "${prefix}.sortmerna.log"
          # Emit empty paired outputs to keep channel arity consistent
          gzip -cn /dev/null > "${prefix}_R1.sortmerna.fastq.gz"
          gzip -cn /dev/null > "${prefix}_R2.sortmerna.fastq.gz"
          exit 0
        fi

        sortmerna \\
            --ref $refs \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir \$TMPDIR/ \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2 \\
            $args

        mv non_rRNA_reads_fwd.f*q.gz "${prefix}_R1.sortmerna.fastq.gz"
        mv non_rRNA_reads_rev.f*q.gz "${prefix}_R2.sortmerna.fastq.gz"
        mv rRNA_reads.log            "${prefix}.sortmerna.log"
        """
    }
}
