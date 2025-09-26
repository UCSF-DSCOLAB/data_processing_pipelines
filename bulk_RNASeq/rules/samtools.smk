# Extract mapped reads (genome/transcriptome), sort/index, stats, CRAM
rule extract_mapped_genome:
    input:
        lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.Aligned.sortedByCoord.out.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.mapped.bam"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools view -@ 1 -b -F 4 -o {output} {input}"

rule extract_mapped_transcriptome:
    input:
        lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.Aligned.toTranscriptome.out.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.transcriptome.mapped.bam"
    shell:
        "samtools view -@ 1 -b -F 4 -o {output} {input}"

rule samtools_sort_mapped:
    input:
        f"{RESULTS_DIR}/star/{{sample}}.mapped.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.bam"
    shell:
        "samtools sort -@ 1 -o {output} -T {wildcards.sample} {input}"

rule samtools_index_bam:
    input:
        f"{RESULTS_DIR}/star/{{sample}}.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.bai"
    shell:
        "samtools index -@ 1 {input} {output}"

rule samtools_stats:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bam",
        bai=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bai"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.stats"
    shell:
        "samtools stats --threads 1 --reference {GENOME} {input.bam} > {output}"

rule samtools_flagstat:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bam",
        bai=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bai"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.flagstat"
    shell:
        "samtools flagstat --threads 1 {input.bam} > {output}"

rule samtools_idxstats:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bam",
        bai=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bai"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.idxstats"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule bam_to_cram_genome:
    input:
        f"{RESULTS_DIR}/star/{{sample}}.mapped.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.mapped.cram"
    shell:
        "samtools view -@ 1 -C --no-PG -T {GENOME} -o {output} {input}"

rule bam_to_cram_transcriptome:
    input:
        f"{RESULTS_DIR}/star/{{sample}}.transcriptome.mapped.bam"
    output:
        f"{RESULTS_DIR}/star/{{sample}}.transcriptome.mapped.cram"
    shell:
        "samtools view -@ 1 -C --no-PG -T {TRANSCRIPT_FASTA} -o {output} {input}"