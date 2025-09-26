# Concatenate per-sample runs (SE/PE)
rule merge_fastq_se:
    input:
        lambda wc: SAMPLES_MAP[wc.sample]["r1"]
    output:
        f"{RESULTS_DIR}/merged_fastq/{{sample}}.merged.fastq.gz"
    shell:
        r"""
        mkdir -p $(dirname {output})
        if [ $(echo {input} | wc -w) -gt 1 ]; then
            cat {input} > {output}
        else
            ln -sf $(realpath {input}) {output}
        fi
        """

rule merge_fastq_pe:
    input:
        r1=lambda wc: SAMPLES_MAP[wc.sample]["r1"],
        r2=lambda wc: SAMPLES_MAP[wc.sample]["r2"]
    output:
        f"{RESULTS_DIR}/merged_fastq/{{sample}}_R1.merged.fastq.gz",
        f"{RESULTS_DIR}/merged_fastq/{{sample}}_R2.merged.fastq.gz"
    shell:
        r"""
        mkdir -p $(dirname {output[0]})
        if [ $(echo {input.r1} | wc -w) -gt 1 ]; then
            cat {input.r1} > {output[0]}
            cat {input.r2} > {output[1]}
        else
            ln -sf $(realpath {input.r1}) {output[0]}
            ln -sf $(realpath {input.r2}) {output[1]}
        fi
        """