process UMITOOLS_PREPARE_FOR_SALMON {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    cpus 2
    memory '31 GB'    
    conda '/c4/home/alaa/miniconda3/envs/umi'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prepare_for_salmon.py \\
        --stdin=$bam \\
        --stdout=${prefix}_presalmon.bam \\
        --log=${prefix}.prepare_for_rsem.log \\
        $args
    """
}