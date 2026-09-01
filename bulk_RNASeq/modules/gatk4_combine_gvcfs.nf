process GATK4_COMBINEGVCFS {
    tag "cohort"
    label 'gatk4_combinegvcfs'
    publishDir "${params.results_directory}/snps/joint", mode: 'copy'
    containerOptions "-B /scratch/"

    input:
    val meta
    path gvcfs
    path tbis
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi"), emit: gvcf

    script:
    def variants = gvcfs.collect { "--variant $it" }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CombineGVCFs \\
        --reference $fasta \\
        $variants \\
        --output cohort.g.vcf.gz \\
        --tmp-dir \$TMPDIR
    """
}
