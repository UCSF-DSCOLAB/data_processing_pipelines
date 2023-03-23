process RSEM_PREPAREREFERENCE {
    publishDir "${params.reference_directory}/genome_dir", mode: 'copy'
    tag "$fasta"
    cpus 16
    memory '64 GB'

    input:
    path fasta, stageAs: "rsem/*"
    path gtf

    output:
    path "rsem"           , emit: index
    path "*transcripts.fa", emit: transcript_fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rsem-prepare-reference \\
        --gtf $gtf \\
        --num-threads $task.cpus \\
        $args \\
        $fasta \\
        rsem/rsem_genome

    mv rsem/rsem_genome.transcripts.fa .
    """
}