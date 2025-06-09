process KALLISTO_INDEX {
    publishDir "${params.reference_directory}", mode: 'copy'
    tag "$fasta"
    label 'kallisto_index'
    memory '64 GB'
    // conda "$baseDir/envs/kallisto.yml"

    input:
    path transcriptome

    output:
    path "*kallisto_index", emit: transcript_index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto index -i kallisto_index $transcriptome
    """
}