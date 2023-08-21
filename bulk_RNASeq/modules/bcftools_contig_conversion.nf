process BCFTOOLS_CONTIG_CONVERSION {
    tag "$meta.id"
    label 'bcftools_contig_conversion'
    cpus 2
    memory '31 GB'
    publishDir "${params.results_directory}/snps", mode: 'copy'
    conda "$baseDir/envs/bcftools.yml"

    input:
    tuple val(meta), path(vcf)
    path format_map

    output:
    tuple val(meta), path("*.formatted.vcf.gz"), emit: formatted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools annotate \\
            --rename-chrs $format_map \\
            --threads $task.cpus -Oz \\
            --output ${prefix}.formatted.vcf.gz \\
            $vcf
            $args
    """
}