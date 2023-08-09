process UMITOOLS_DEDUP {
    tag "$meta.id"
    cpus 2
    memory '31 GB'
    conda '/c4/home/alaa/miniconda3/envs/umi'

    input:
    tuple val(meta), path(bam), path(bai)
    val get_output_stats

    output:
    tuple val(meta), path("*.bam")             , emit: bam
    tuple val(meta), path("*edit_distance.tsv"), optional:true, emit: tsv_edit_distance
    tuple val(meta), path("*per_umi.tsv")      , optional:true, emit: tsv_per_umi
    tuple val(meta), path("*per_position.tsv") , optional:true, emit: tsv_umi_per_position

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    def stats = get_output_stats ? "--output-stats $prefix" : ""

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        --stdin $bam \\
        --stdout ${prefix}.bam \\
        $stats \\
        $paired \\
        $args
    """
}