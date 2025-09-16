# Extract mapped reads (genome/transcriptome), sort/index, stats, CRAM
rule extract_mapped_genome:
    input:
        lambda wc: star_sorted_bam(wc.sample)
    output:
        lambda wc: mapped_bam(wc.sample)
    conda:
        "envs/samtools.yml"
    shell:
        "samtools view -@ 1 -b -F 4 -o {output} {input}"

rule extract_mapped_transcriptome:
    input:
        lambda wc: star_transcriptome_bam(wc.sample)
    output:
        lambda wc: mapped_transcriptome_bam(wc.sample)
    conda:
        "envs/samtools.yml"
    shell:
        "samtools view -@ 1 -b -F 4 -o {output} {input}"

rule samtools_sort_mapped:
    input:
        lambda wc: mapped_bam(wc.sample)
    output:
        lambda wc: sorted_bam(wc.sample)
    conda:
        "envs/samtools.yml"
    shell:
        "samtools sort -@ 1 -o {output} -T {wildcards.sample} {input}"

rule samtools_index_bam:
    input:
        lambda wc: sorted_bam(wc.sample)
    output:
        lambda wc: sorted_bai(wc.sample)
    conda:
        "envs/samtools.yml"
    shell:
        "samtools index -@ 1 {input}"

rule samtools_stats:
    input:
        bam=lambda wc: sorted_bam(wc.sample),
        bai=lambda wc: sorted_bai(wc.sample)
    output:
        f"{RESULTS_DIR}/star/{{sample}}.stats"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools stats --threads 1 --reference {GENOME} {input.bam} > {output}"

rule samtools_flagstat:
    input:
        bam=lambda wc: sorted_bam(wc.sample),
        bai=lambda wc: sorted_bai(wc.sample)
    output:
        f"{RESULTS_DIR}/star/{{sample}}.flagstat"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools flagstat --threads 1 {input.bam} > {output}"

rule samtools_idxstats:
    input:
        bam=lambda wc: sorted_bam(wc.sample),
        bai=lambda wc: sorted_bai(wc.sample)
    output:
        f"{RESULTS_DIR}/star/{{sample}}.idxstats"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule bam_to_cram_genome:
    input:
        lambda wc: mapped_bam(wc.sample)
    output:
        f"{RESULTS_DIR}/star/{{sample}}.mapped.cram"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools view -@ 1 -C --no-PG -T {GENOME} -o {output} {input}"

rule bam_to_cram_transcriptome:
    input:
        lambda wc: mapped_transcriptome_bam(wc.sample)
    output:
        f"{RESULTS_DIR}/star/{{sample}}.transcriptome.mapped.cram"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools view -@ 1 -C --no-PG -T {TRANSCRIPT_FASTA} -o {output} {input}"