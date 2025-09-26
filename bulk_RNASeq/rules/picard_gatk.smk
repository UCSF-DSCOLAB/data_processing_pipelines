# MarkDuplicates, SplitNCigarReads, BQSR, ApplyBQSR, HaplotypeCaller, VariantFiltration
rule picard_markduplicates:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bam",
        bai=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.bai"
    output:
        bam=f"{RESULTS_DIR}/star/{{sample}}.picard.bam",
        bai=f"{RESULTS_DIR}/star/{{sample}}.picard.bai",
        metrics=f"{RESULTS_DIR}/star/{{sample}}.MarkDuplicates.metrics.txt"
    shell:
        r"""
        gatk --java-options "-Xmx4g" MarkDuplicates \
          --INPUT {input.bam} \
          --OUTPUT {output.bam} \
          --REFERENCE_SEQUENCE {GENOME} \
          --METRICS_FILE {output.metrics}
        samtools index {output.bam} {output.bai}
        """

rule gatk_splitncigarreads:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.picard.bam" if False else f"{RESULTS_DIR}/star/{wc.sample}.picard.bam",
        bai=lambda wc: f"{RESULTS_DIR}/star/{wc.sample}.picard.bai"
    output:
        bam=f"{RESULTS_DIR}/snps/{{sample}}.split.bam",
        bai=f"{RESULTS_DIR}/snps/{{sample}}.split.bam.bai"
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
        bam=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.split.bam",
        bai=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.split.bam.bai"
    output:
        f"{RESULTS_DIR}/snps/{{sample}}.table"
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
        bam=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.split.bam",
        bai=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.split.bam.bai",
        table=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.table"
    output:
        bam=f"{RESULTS_DIR}/snps/{{sample}}_bqsr.bam"
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
        f"{RESULTS_DIR}/snps/{{sample}}_bqsr.bam"
    output:
        f"{RESULTS_DIR}/snps/{{sample}}_bqsr.bam.bai"
    shell:
        "samtools index {input} {output}"

rule gatk_haplotypecaller:
    input:
        bam=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}_bqsr.bam",
        bai=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}_bqsr.bam.bai"
    output:
        vcf=f"{RESULTS_DIR}/snps/{{sample}}.vcf.gz",
        tbi=f"{RESULTS_DIR}/snps/{{sample}}.vcf.gz.tbi"
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
        vcf=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.vcf.gz",
        tbi=lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.vcf.gz.tbi"
    output:
        vcf=f"{RESULTS_DIR}/snps/{{sample}}.filtered.vcf.gz",
        tbi=f"{RESULTS_DIR}/snps/{{sample}}.filtered.vcf.gz.tbi"
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