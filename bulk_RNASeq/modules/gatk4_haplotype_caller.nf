process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    publishDir "${params.results_directory}/snps", mode: 'copy'

    input:
    tuple val(meta), path(input), path(input_index)
    path  fasta
    path  fai
    path  dict
    path  known_sites
    path  known_sites_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , optional:true, emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_command = "--reference $fasta"
    def dbsnp_command = known_sites ? "--dbsnp $known_sites" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        --input $input \\
        --output ${prefix}.vcf.gz \\
        $reference_command \\
        $dbsnp_command \\
        --tmp-dir $params.tmp_dir \\
        $args
    """
}
