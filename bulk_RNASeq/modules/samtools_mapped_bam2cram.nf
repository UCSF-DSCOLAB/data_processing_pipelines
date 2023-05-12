process EXTRACT_MAPPED_READS_BAM2CRAM {
    tag "$meta.id"
    cpus 32
    memory '64 GB'
    publishDir "${params.results_directory}/star", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.trimmed.star.mapped.bam'), emit: mapped_bam
    path('*.trimmed.star.mapped.cram')


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}."
    """
    samtools view  -b \\
        -F 0x4 \\
        -@ ${task.cpus} \\
        --no-PG \\
        -o ${prefix}.trimmed.star.Transcriptome.mapped.bam \\
        ${bam}

    samtools view -@ ${task.cpus} \\
        -C \\
        --no-PG \\
        -T ${params.ref.rsem_star_dir}/${params.ref.genome_version}.transcripts.fa \\
        -o ${prefix}.trimmed.star.mapped.cram \\
        ${prefix}.trimmed.star.mapped.bam
    """
}