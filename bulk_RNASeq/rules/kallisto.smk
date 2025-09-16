# Kallisto quant and merged counts
rule kallisto:
    input:
        reads=lambda wc: reads_for_downstream(wc.sample),
        index=KALLISTO_INDEX
    output:
        f"{RESULTS_DIR}/kallisto/{{sample}}.kallisto.abundance.tsv",
        f"{RESULTS_DIR}/kallisto/{{sample}}.kallisto.abundance.h5",
        f"{RESULTS_DIR}/kallisto/{{sample}}.kallisto.run_info.json",
        f"{RESULTS_DIR}/kallisto/{{sample}}.kallisto.log"
    threads: 4
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/kallisto
        READS=({input.reads})
        if [ ${{#READS[@]}} -eq 1 ]; then
            (kallisto quant --index {input.index} --output-dir {RESULTS_DIR}/kallisto/{wildcards.sample} \
                -t {threads} --single -l {FRAG_MEAN} -s {FRAG_STD} ${READS[0]} \
                > {RESULTS_DIR}/kallisto/{wildcards.sample}.kallisto.log)
        else
            (kallisto quant --index {input.index} --output-dir {RESULTS_DIR}/kallisto/{wildcards.sample} \
                -t {threads} ${READS[0]} ${READS[1]} \
                > {RESULTS_DIR}/kallisto/{wildcards.sample}.kallisto.log)
        fi
        mv {RESULTS_DIR}/kallisto/{wildcards.sample}/abundance.h5 {RESULTS_DIR}/kallisto/{wildcards.sample}.kallisto.abundance.h5
        mv {RESULTS_DIR}/kallisto/{wildcards.sample}/abundance.tsv {RESULTS_DIR}/kallisto/{wildcards.sample}.kallisto.abundance.tsv
        mv {RESULTS_DIR}/kallisto/{wildcards.sample}/run_info.json {RESULTS_DIR}/kallisto/{wildcards.sample}.kallisto.run_info.json
        rmdir {RESULTS_DIR}/kallisto/{wildcards.sample}
        """

rule merge_kallisto_counts:
    input:
        lambda wc: [f"{RESULTS_DIR}/kallisto/{s}.kallisto.abundance.tsv" for s in SAMPLES]
    output:
        f"{RESULTS_DIR}/merged_results/merged_counts.tsv"
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/merged_results
        merge_kallisto_counts.py "{','.join(input)}" {output}
        """