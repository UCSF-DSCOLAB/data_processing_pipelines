# Optional contig format conversion, sort, index, merge
rule bcftools_contig_conversion:
    input:
        lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.filtered.vcf.gz"
    output:
        f"{RESULTS_DIR}/snps/{{sample}}.formatted.vcf.gz"
    params:
        mapfile=CONTIG_FORMAT_MAP
    run:
        if not FORMAT_CONTIGS or not CONTIG_FORMAT_MAP:
            shell("ln -sf {input} {output}")
        else:
            shell("bcftools annotate --rename-chrs {params.mapfile} -Oz -o {output} {input}")

rule bcftools_sort_vcf:
    input:
        lambda wc: f"{RESULTS_DIR}/snps/{wc.sample}.formatted.vcf.gz" if FORMAT_CONTIGS and CONTIG_FORMAT_MAP else f"{RESULTS_DIR}/snps/{wc.sample}.filtered.vcf.gz"
    output:
        f"{RESULTS_DIR}/snps/{{sample}}.sorted.vcf.gz"
    shell:
        "bcftools sort -Oz -o {output} {input}"

rule bcftools_index_vcf:
    input:
        f"{RESULTS_DIR}/snps/{{sample}}.sorted.vcf.gz"
    output:
        f"{RESULTS_DIR}/snps/{{sample}}.sorted.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"

rule bcftools_merge:
    input:
        vcfs=lambda wc: [f"{RESULTS_DIR}/snps/{s}.sorted.vcf.gz" for s in SAMPLES]
    output:
        f"{RESULTS_DIR}/merged_results/merged_snps.bcf"
    threads: 8
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/merged_results
        if [ $(echo {input.vcfs} | wc -w) -gt 1 ]; then
            bcftools merge --threads {threads} --output-type b --output {output} {input.vcfs}
        else
            bcftools view --output-type b --output {output} {input.vcfs}
        fi
        """