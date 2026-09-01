process GATK4_SELECT_SNPS {
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
    tuple val(meta), path("cohort.snps.vcf.gz"), path("cohort.snps.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SelectVariants \\
        --reference $fasta \\
        --variant $vcf \\
        --select-type-to-include SNP \\
        --restrict-alleles-to BIALLELIC \\
        --output cohort.snps.vcf.gz \\
        --tmp-dir \$TMPDIR
    """
}
