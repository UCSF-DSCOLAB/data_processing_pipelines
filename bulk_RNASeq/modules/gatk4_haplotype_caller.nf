process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'gatk4_haplotypecaller', 'per_sample'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    memory {
        // File size in GB
        fileSize = input.size() / (1024 * 1024 * 1024)
        return 17.GB + (1.GB * fileSize * 3)
    }

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
    def soft_clipped = params.gatk_dont_use_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    def min_conf = "--standard-min-confidence-threshold-for-calling ${params.gatk_standard_min_confidence}"
    def min_pruning = "--min-pruning ${params.gatk_min_pruning}"
    def recover_branches = params.gatk_recover_all_dangling_branches ? "--recover-all-dangling-branches true" : ""
    def allow_nonunique = params.gatk_allow_nonunique_kmer ? "--allow-nonuniquekmer true" : ""
    def max_mnp_distance = "--max-mnp-distance ${params.gatk_max_mnp_distance}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        --input $input \\
        --output ${prefix}.vcf.gz \\
        $reference_command \\
        $dbsnp_command \\
        --tmp-dir \$PWD \\
        $soft_clipped \\
        $min_conf \\
        $min_pruning \\
        $recover_branches \\
        $allow_nonunique \\
        $max_mnp_distance \\
        $args
    """
}
