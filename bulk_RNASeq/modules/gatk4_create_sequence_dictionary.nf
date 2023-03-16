process GATK4_GENERATE_SEQUENCE_DICTIONARY {
    publishDir "${params.reference_directory}/", mode: 'copy'
    cpus 2
    memory '314 GB'

    input:
    path  genome

    output:
    path  "*.dict", emit: genome_dict

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CreateSequenceDictionary \\
         --REFERENCE $genome \\
         --URI $genome \\
         --TMP_DIR "${params.tmp_dir}/gatk4_create_sequence_dictionary"
        $args
    """
}