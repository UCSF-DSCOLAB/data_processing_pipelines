process STAR_INDEX_GENOME {
    publishDir "${params.reference_directory}/genome_dir", mode: 'copy'
    cpus 92
    memory '314 GB'

    input:
    path  genome
    path  gtf

    output:
    path "genome_dir"         , emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def memory   = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def args = task.ext.args ?: ''
    """
    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles $genome \\
         --sjdbGTFfile $gtf \\
         $memory \\
         --runThreadN $task.cpus \\
         --outTmpDir "${params.tmp_dir}/star_index_genome"
         $args
    """
}