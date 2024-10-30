process SALMON_QUANT {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    publishDir "${params.results_directory}/salmon", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path  gtf
    path  transcript_fasta
    val   lib_type

    output:
    tuple val(meta), path("${prefix}") , emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def reference   = "-t $transcript_fasta"
    def input_reads = "-a $bam"

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $args \\
        -o $prefix

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi
    """
}