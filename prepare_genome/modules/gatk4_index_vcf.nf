process GATK4_INDEX_VCF {
    label "gatk_index_vcf"
    publishDir "${params.reference_directory}/", mode: 'copy'
    memory '256 GB'
    // conda "$baseDir/envs/gatk.yml"

    input:
    path  vcf

    output:
    path  "*.tbi", emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk --java-options "-Xmx314g" IndexFeatureFile \\
    -I $vcf
    """
}