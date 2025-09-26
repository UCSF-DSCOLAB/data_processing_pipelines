# Validation rule (optional materialization)
rule validate_samplesheet:
    input:
        lambda wc: SAMPLESHEET
    output:
        f"{RESULTS_DIR}/samplesheet.valid.csv"
    shell:
        "validate_sample_sheet.py {input} {output}"