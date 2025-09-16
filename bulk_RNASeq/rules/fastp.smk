# Adapter and quality trimming
rule fastp:
    input:
        lambda wc: merged_fastq_se(wc.sample) if SAMPLES_MAP[wc.sample]["se"] else merged_fastq_pe(wc.sample)
    output:
        lambda wc: fastp_outputs_se(wc.sample) if SAMPLES_MAP[wc.sample]["se"] else fastp_outputs_pe(wc.sample)
    threads: 6
    conda:
        "envs/fastp.yml"
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/trimmed_reads
        if [ {1 if SAMPLES_MAP[wildcards.sample]["se"] else 0} -eq 1 ]; then
          fastp \
            --in1 {input} \
            --out1 {output[0]} \
            --length_required 20 \
            {("--adapter_sequence " + ADAPTER_1) if ADAPTER_1 else ""} \
            {("--adapter_sequence_r2 " + ADAPTER_2) if ADAPTER_2 else ""} \
            --correction \
            --trim_poly_g \
            --thread {threads} \
            --json {output[1]} \
            --html {output[2]}
        else
          fastp \
            --in1 {input[0]} \
            --in2 {input[1]} \
            --out1 {output[0]} \
            --out2 {output[1]} \
            --length_required 20 \
            {("--adapter_sequence " + ADAPTER_1) if ADAPTER_1 else ""} \
            {("--adapter_sequence_r2 " + ADAPTER_2) if ADAPTER_2 else ""} \
            --correction \
            --trim_poly_g \
            --thread {threads} \
            --json {output[2]} \
            --html {output[3]}
        fi
        """