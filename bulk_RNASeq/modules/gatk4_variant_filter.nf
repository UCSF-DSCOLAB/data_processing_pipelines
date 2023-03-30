process GATK4_VARIANTFILTRATION {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    publishDir "${params.results_directory}/snps", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration \\
        --variant $vcf \\
        --cluster $params.gatk_vf_cluster_size \\
        --filter-name FS -filter "FS > $params.gatk_vf_fs_filter" \\
        --filter-name QD -filter "QD < $params.gatk_vf_qd_filter" \\
        --reference $fasta \\
        --window $params.gatk_vf_window_size \\
        --output ${prefix}.filtered.vcf.gz \\
        --tmp-dir . \\
        $args
    """
}