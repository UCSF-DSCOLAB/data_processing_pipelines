process GATK4_GENOTYPEGVCFS {
    tag "cohort"
    label 'gatk4_genotypegvcfs'
    publishDir "${params.results_directory}/snps/joint", mode: 'copy'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("cohort.raw.vcf.gz"), path("cohort.raw.vcf.gz.tbi"), emit: vcf

    script:
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        --reference $fasta \\
        --variant $gvcf \\
        $dbsnp_command \\
        --output cohort.raw.vcf.gz \\
        --tmp-dir \$TMPDIR
    """
}
