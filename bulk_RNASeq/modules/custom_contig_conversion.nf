process CUSTOM_CONTIG_CONVERSION {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    conda '/c4/home/alaa/miniconda3/envs/umi'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.formatted.vcf.gz"), emit: formatted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools annotate \\
            --rename-chrs $params.contig_format_map \\
            --threads $task.cpus -Oz \\
            --output ${prefix}.formatted.vcf.gz \\
            $vcf
            $args
    """
}