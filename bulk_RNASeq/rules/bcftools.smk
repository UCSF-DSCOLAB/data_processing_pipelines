# Optional contig format conversion, sort, index, merge
rule bcftools_contig_conversion:
    input:
        lambda wc: filtered_vcf(wc.sample)
    output:
        lambda wc: formatted_vcf(wc.sample)
    conda:
        "envs/bcftools.yml"
    params:
        mapfile=CONTIG_FORMAT_MAP
    run:
        if not FORMAT_CONTIGS or not CONTIG_FORMAT_MAP:
            shell("ln -sf {input} {output}")
        else:
            shell("bcftools annotate --rename-chrs {params.mapfile} -Oz -o {output} {input}")

rule bcftools_sort_vcf:
    input:
        lambda wc: formatted_vcf(wc.sample) if FORMAT_CONTIGS and CONTIG_FORMAT_MAP else filtered_vcf(wc.sample)
    output:
        lambda wc: sorted_vcf(wc.sample)
    conda:
        "envs/bcftools.yml"
    shell:
        "bcftools sort -Oz -o {output} {input}"

rule bcftools_index_vcf:
    input:
        lambda wc: sorted_vcf(wc.sample)
    output:
        lambda wc: sorted_tbi(wc.sample)
    conda:
        "envs/bcftools.yml"
    shell:
        "bcftools index --tbi {input}"

rule bcftools_merge:
    input:
        vcfs=lambda wc: [sorted_vcf(s) for s in SAMPLES]
    output:
        f"{RESULTS_DIR}/merged_results/merged_snps.bcf"
    threads: 8
    conda:
        "envs/bcftools.yml"
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/merged_results
        if [ $(echo {input.vcfs} | wc -w) -gt 1 ]; then
            bcftools merge --threads {threads} --output-type b --output {output} {input.vcfs}
        else
            bcftools view --output-type b --output {output} {input.vcfs}
        fi
        """