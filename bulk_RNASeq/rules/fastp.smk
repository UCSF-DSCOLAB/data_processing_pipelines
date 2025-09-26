# Adapter and quality trimming (split into SE/PE to keep outputs static)

rule fastp_se:
    input:
        lambda wc: merged_fastq_se(wc.sample)
    output:
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}.fastp.fastq.gz",
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}.fastp.json",
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}.fastp.html"
    threads: 6
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/trimmed_reads
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
        """

rule fastp_pe:
    input:
        lambda wc: merged_fastq_pe(wc.sample)
    output:
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R1.fastp.fastq.gz",
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R2.fastp.fastq.gz",
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}.fastp.json",
        f"{RESULTS_DIR}/trimmed_reads/{{sample}}.fastp.html"
    threads: 6
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/trimmed_reads
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
        """