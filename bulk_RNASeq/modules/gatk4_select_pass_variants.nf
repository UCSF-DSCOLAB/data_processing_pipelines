process GATK4_SELECT_PASS_VARIANTS {
    tag "cohort"
    label 'gatk4_selectvariants'
    publishDir "${params.results_directory}/snps/filtered", mode: 'copy'
    containerOptions "-B /scratch/"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("cohort.pass.vcf.gz"), path("cohort.pass.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SelectVariants \\
        --reference $fasta \\
        --variant $vcf \\
        --exclude-filtered true \\
        --output cohort.pass.vcf.gz \\
        --tmp-dir \$TMPDIR
    """
}
