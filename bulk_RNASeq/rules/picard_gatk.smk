# MarkDuplicates, SplitNCigarReads, BQSR, ApplyBQSR, HaplotypeCaller, VariantFiltration
rule picard_markduplicates:
    input:
        bam=lambda wc: sorted_bam(wc.sample),
        bai=lambda wc: sorted_bai(wc.sample)
    output:
        bam=lambda wc: picard_bam(wc.sample),
        bai=lambda wc: picard_bai(wc.sample),
        metrics=lambda wc: picard_metrics(wc.sample)
    conda:
        "envs/gatk.yml"
    shell:
        r"""
        gatk --java-options "-Xmx4g" MarkDuplicates \
          --INPUT {input.bam} \
          --OUTPUT {output.bam} \
          --REFERENCE_SEQUENCE {GENOME} \
          --METRICS_FILE {output.metrics}
        samtools index {output.bam}
        """

rule gatk_splitncigarreads:
    input:
        bam=lambda wc: picard_bam(wc.sample),
        bai=lambda wc: picard_bai(wc.sample)
    output:
        bam=lambda wc: split_bam(wc.sample),
        bai=lambda wc: split_bai(wc.sample)
    conda:
        "envs/gatk.yml"
    shell:
        r"""
        gatk --java-options "-Xmx16g" SplitNCigarReads \
          --input {input.bam} \
          --output {output.bam} \
          --reference {GENOME}
        samtools index {output.bam}
        """

rule gatk_base_recalibrator:
    input:
        bam=lambda wc: split_bam(wc.sample),
        bai=lambda wc: split_bai(wc.sample)
    output:
        lambda wc: bqsr_table(wc.sample)
    conda:
        "envs/gatk.yml"
    params:
        known_sites=(f" --known-sites {DBSNP}" if DBSNP else "")
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/snps
        gatk --java-options "-Xmx8g" BaseRecalibrator \
          --input {input.bam} \
          --output {output} \
          --reference {GENOME} {params.known_sites}
        """

rule gatk_apply_bqsr:
    input:
        bam=lambda wc: split_bam(wc.sample),
        bai=lambda wc: split_bai(wc.sample),
        table=lambda wc: bqsr_table(wc.sample)
    output:
        bam=lambda wc: bqsr_bam(wc.sample)
    conda:
        "envs/gatk.yml"
    shell:
        r"""
        gatk --java-options "-Xmx8g" ApplyBQSR \
          --input {input.bam} \
          --output {output.bam} \
          --reference {GENOME} \
          --bqsr-recal-file {input.table}
        """

rule samtools_index_bqsr:
    input:
        lambda wc: bqsr_bam(wc.sample)
    output:
        lambda wc: bqsr_bai(wc.sample)
    conda:
        "envs/samtools.yml"
    shell:
        "samtools index {input}"

rule gatk_haplotypecaller:
    input:
        bam=lambda wc: bqsr_bam(wc.sample),
        bai=lambda wc: bqsr_bai(wc.sample)
    output:
        vcf=lambda wc: haplotype_vcf(wc.sample),
        tbi=lambda wc: haplotype_tbi(wc.sample)
    conda:
        "envs/gatk.yml"
    params:
        soft_clip=("--dont-use-soft-clipped-bases true" if GATK_PARAMS["dont_use_soft_clipped_bases"] else ""),
        min_conf=f"--standard-min-confidence-threshold-for-calling {GATK_PARAMS['standard_min_conf']}",
        min_pruning=f"--min-pruning {GATK_PARAMS['min_pruning']}",
        recover=("--recover-all-dangling-branches true" if GATK_PARAMS["recover_all_dangling_branches"] else ""),
        allow_nonuniq=("--allow-nonuniquekmer true" if GATK_PARAMS["allow_nonuniquekmer"] else ""),
        max_mnp=f"--max-mnp-distance {GATK_PARAMS['max_mnp_distance']}",
        dbsnp=(f"--dbsnp {DBSNP}" if DBSNP else "")
    shell:
        r"""
        gatk --java-options "-Xmx16g" HaplotypeCaller \
          --input {input.bam} \
          --output {output.vcf} \
          --reference {GENOME} \
          {params.dbsnp} \
          {params.soft_clip} \
          {params.min_conf} \
          {params.min_pruning} \
          {params.recover} \
          {params.allow_nonuniq} \
          {params.max_mnp}
        """

rule gatk_variantfiltration:
    input:
        vcf=lambda wc: haplotype_vcf(wc.sample),
        tbi=lambda wc: haplotype_tbi(wc.sample)
    output:
        vcf=lambda wc: filtered_vcf(wc.sample),
        tbi=lambda wc: filtered_tbi(wc.sample)
    conda:
        "envs/gatk.yml"
    shell:
        r"""
        gatk --java-options "-Xmx8g" VariantFiltration \
          --variant {input.vcf} \
          --cluster {GATK_PARAMS['vf_cluster_size']} \
          --filter-name FS -filter "FS > {GATK_PARAMS['vf_fs_filter']}" \
          --filter-name QD -filter "QD < {GATK_PARAMS['vf_qd_filter']}" \
          --reference {GENOME} \
          --window {GATK_PARAMS['vf_window_size']} \
          --output {output.vcf}
        """