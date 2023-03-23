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

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , optional:true, emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        --input $input \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        --tmp-dir $params.tmp_dir \\
        $args
    """
}
