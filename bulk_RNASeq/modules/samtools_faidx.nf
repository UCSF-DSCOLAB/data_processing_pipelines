process SAMTOOLS_FAIDX {
    publishDir "${params.reference_directory}/", mode: 'copy'
    tag "$genome"
    cpus 2
    memory '314 GB'

    input:
    path(genome)

    output:
    path ("*.fai"), emit: genome_samtools_index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools \\
        faidx \\
        $genome
    """
}