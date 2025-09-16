# Final targets and DAG connections
MERGED_COUNTS = f"{RESULTS_DIR}/merged_results/merged_counts.tsv"
MERGED_BCF    = f"{RESULTS_DIR}/merged_results/merged_snps.bcf"
MULTIQC_HTML  = f"{RESULTS_DIR}/multiqc/multiqc_report.html"

rule all:
    input:
        # Trigger SE/PE concatenation
        [merged_fastq_se(s) for s in SAMPLES_SE] +
        [p for s in SAMPLES_PE for p in merged_fastq_pe(s)] +
        # fastp outputs
        [fastp_outputs_se(s)[0] for s in SAMPLES_SE] +
        [p for s in SAMPLES_PE for p in fastp_outputs_pe(s)[0:2]] +
        # Optional SortMeRNA
        ([] if not USE_RRNA_FILTER else
         [f"{RESULTS_DIR}/trimmed_cleaned_reads/{s}.sortmerna.fastq.gz" for s in SAMPLES_SE] +
         [f"{RESULTS_DIR}/trimmed_cleaned_reads/{s}_R1.sortmerna.fastq.gz" for s in SAMPLES_PE] +
         [f"{RESULTS_DIR}/trimmed_cleaned_reads/{s}_R2.sortmerna.fastq.gz" for s in SAMPLES_PE]) +
        # STAR outputs
        [star_sorted_bam(s) for s in SAMPLES] +
        [star_transcriptome_bam(s) for s in SAMPLES] +
        [star_gene_counts(s) for s in SAMPLES] +
        [star_log_final(s) for s in SAMPLES] +
        # Mapped extraction, sorting, indexing, stats
        [mapped_bam(s) for s in SAMPLES] +
        [sorted_bam(s) for s in SAMPLES] +
        [sorted_bai(s) for s in SAMPLES] +
        [f"{RESULTS_DIR}/star/{s}.stats" for s in SAMPLES] +
        [f"{RESULTS_DIR}/star/{s}.flagstat" for s in SAMPLES] +
        [f"{RESULTS_DIR}/star/{s}.idxstats" for s in SAMPLES] +
        # CRAMs
        [f"{RESULTS_DIR}/star/{s}.mapped.cram" for s in SAMPLES] +
        [f"{RESULTS_DIR}/star/{s}.transcriptome.mapped.cram" for s in SAMPLES] +
        # Picard/GATK
        [picard_bam(s) for s in SAMPLES] +
        [split_bam(s) for s in SAMPLES] +
        [bqsr_table(s) for s in SAMPLES] +
        [bqsr_bam(s) for s in SAMPLES] +
        [bqsr_bai(s) for s in SAMPLES] +
        # VCFs
        [haplotype_vcf(s) for s in SAMPLES] +
        [filtered_vcf(s) for s in SAMPLES] +
        [sorted_vcf(s) for s in SAMPLES] +
        [sorted_tbi(s) for s in SAMPLES] +
        # Kallisto
        [kallisto_outputs(s)["abundance_tsv"] for s in SAMPLES] +
        [kallisto_outputs(s)["abundance_h5"] for s in SAMPLES] +
        [kallisto_outputs(s)["run_info"] for s in SAMPLES] +
        [kallisto_outputs(s)["log"] for s in SAMPLES] +
        # Merged results
        MERGED_COUNTS,
        MERGED_BCF,
        MULTIQC_HTML