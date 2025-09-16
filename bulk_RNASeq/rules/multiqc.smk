# MultiQC aggregation
rule multiqc:
    input:
        lambda wc: multiqc_inputs()
    output:
        f"{RESULTS_DIR}/multiqc/multiqc_report.html"
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/multiqc
        multiqc -f {RESULTS_DIR} -o {RESULTS_DIR}/multiqc
        """