process SAMTOOLS_FAIDX {
    label 'samtools_faidx'
    publishDir "${params.reference_directory}/", mode: 'copy'
    tag "$genome"
    memory '256 GB'
    // conda "$baseDir/envs/samtools.yml"

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