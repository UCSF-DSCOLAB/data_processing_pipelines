process STAR_INDEX_GENOME {
    label 'star_index_genome'
    publishDir "${params.reference_directory}", mode: 'copy'
    memory '256 GB'
    // conda "$baseDir/envs/star.yml"

    input:
    path  genome
    path  gtf

    output:
    path "star_index"         , emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def memory   = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def args = task.ext.args ?: ''
    """
    rm -rf "${params.tmp_dir}/star_index"
    mkdir "${params.tmp_dir}"

    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles $genome \\
         --sjdbGTFfile $gtf \\
         $memory \\
         --runThreadN $task.cpus \\
         --outTmpDir "${params.tmp_dir}/star_index"
         $args
    """
}