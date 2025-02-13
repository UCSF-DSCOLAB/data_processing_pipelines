process RSEM_PREPAREREFERENCE {
    label 'rsem_prepare_reference'
    publishDir "${params.reference_directory}", mode: 'copy'
    tag "$fasta"
    // cpus 32
    memory '64 GB'
    // conda "$baseDir/envs/rsem.yml"

    input:
    path fasta, stageAs: "rsem_index/*"
    path gtf

    output:
    path "rsem_index"           , emit: index
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
        rsem_index/rsem_genome

    mv rsem_index/rsem_genome.transcripts.fa .
    """
}