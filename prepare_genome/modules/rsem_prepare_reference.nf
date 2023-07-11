process RSEM_PREPAREREFERENCE {
    publishDir "${params.reference_directory}", mode: 'copy'
    tag "$fasta"
    cpus 32
    memory '64 GB'
    conda "$baseDir/envs/rsem.yml"

    input:
    path fasta, stageAs: "rsem/*"
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
        rsem/rsem_genome

    mv rsem/rsem_genome.transcripts.fa .
    """
}