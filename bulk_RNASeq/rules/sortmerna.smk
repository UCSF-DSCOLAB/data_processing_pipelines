# Optional rRNA filtering with SortMeRNA (SE/PE separately)
# Only used downstream if filter_rrna is true.

rule sortmerna_se:
    input:
        lambda wc: fastp_outputs_se(wc.sample)[0]
    output:
        f"{RESULTS_DIR}/trimmed_cleaned_reads/{{sample}}.sortmerna.fastq.gz",
        f"{RESULTS_DIR}/trimmed_cleaned_reads/{{sample}}.sortmerna.log"
    threads: 6
    params:
        refs=rrna_ref_args()
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/trimmed_cleaned_reads
        sortmerna \
            {params.refs} \
            --reads {input} \
            --threads {threads} \
            --workdir {TMP_DIR}/ \
            --aligned rRNA_reads \
            --fastx \
            --other non_rRNA_reads
        mv non_rRNA_reads.f*q.gz {output[0]}
        mv rRNA_reads.log {output[1]}
        """

rule sortmerna_pe:
    input:
        lambda wc: fastp_outputs_pe(wc.sample)[0:2]
    output:
        f"{RESULTS_DIR}/trimmed_cleaned_reads/{{sample}}_R1.sortmerna.fastq.gz",
        f"{RESULTS_DIR}/trimmed_cleaned_reads/{{sample}}_R2.sortmerna.fastq.gz",
        f"{RESULTS_DIR}/trimmed_cleaned_reads/{{sample}}.sortmerna.log"
    threads: 6
    params:
        refs=rrna_ref_args()
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/trimmed_cleaned_reads
        sortmerna \
            {params.refs} \
            --reads {input[0]} \
            --reads {input[1]} \
            --threads {threads} \
            --workdir {TMP_DIR}/ \
            --aligned rRNA_reads \
            --fastx \
            --other non_rRNA_reads \
            --paired_in \
            --out2
        mv non_rRNA_reads_fwd.f*q.gz {output[0]}
        mv non_rRNA_reads_rev.f*q.gz {output[1]}
        mv rRNA_reads.log {output[2]}
        """