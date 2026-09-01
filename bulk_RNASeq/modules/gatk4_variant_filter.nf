process GATK4_VARIANTFILTRATION {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'gatk4_variantfiltration', 'per_sample'
    publishDir "${params.results_directory}/snps/filtered", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcf.size() / (1024 * 1024 * 1024)
        return 5.GB + (1.GB * fileSize)
    }

    containerOptions "-B /scratch/"

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
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantFiltration \\
        --variant $vcf \\
        --cluster $params.gatk_vf_cluster_size \\
        --filter-name FS -filter "FS > $params.gatk_vf_fs_filter" \\
        --filter-name QD -filter "QD < $params.gatk_vf_qd_filter" \\
        --filter-name QUAL -filter "QUAL < $params.gatk_vf_qual_filter" \\
        --filter-name SOR -filter "SOR > $params.gatk_vf_sor_filter" \\
        --filter-name MQ -filter "MQ < $params.gatk_vf_mq_filter" \\
        --filter-name MQRankSum -filter "MQRankSum < $params.gatk_vf_mq_rank_sum_filter" \\
        --filter-name ReadPosRankSum -filter "ReadPosRankSum < $params.gatk_vf_read_pos_rank_sum_filter" \\
        --genotype-filter-name LowDP --genotype-filter-expression "DP < $params.genotype_min_depth" \\
        --genotype-filter-name LowGQ --genotype-filter-expression "GQ < $params.genotype_min_gq" \\
        --set-filtered-genotype-to-no-call true \\
        --reference $fasta \\
        --window $params.gatk_vf_window_size \\
        --output ${prefix}.filtered.vcf.gz \\
        --tmp-dir \$TMPDIR \\
        $args
    """
}
