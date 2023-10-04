process BCFTOOLS_MERGE_VCF {
    tag "$meta.id"
    label 'bcftools_merge_vcf'
    publishDir "${params.results_directory}/merged_results", mode: 'copy'
    memory {
        // File size in GB
        fileSize = vcfs.size() / (1024 * 1024 * 1024)
        return 5.GB * (1 + (fileSize * 0.01))
    }

    input:
    val meta 
    path vcfs
    path tbis

    output:
    path("*.{bcf,vcf}{,.gz}"), emit: merged_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create a comma-separated list of VCF file names
    bcftools merge \\
    --threads $task.cpus \\
    --output-type b \\
    --output merged_snps.bcf \\
    ${vcfs.join(' ')}
    """
}