process STAR_ALIGN {
    tag "$meta.id"
    cpus 2
    memory '31 GB'

    input:
    tuple val(meta), path(reads)
    path gtf
    path genomeDir
    path tmpDir

    output:
    tuple val(meta), path('*d.sortedByCoord.out.bam'), emit: bam
    tuple val(meta), path('*d.sortedByCoord.out.bam.bai'), emit: bai
    tuple val(meta), path('*d.toTranscriptome.out.bam'), emit: transcriptome_bam


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}."
    """
    STAR \\
        --readFilesIn $reads  \\
        --genomeDir $genomeDir \\
        --runThreadN $task.cpus \\
        --sjdbGTFfile $gtf \\
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --quantTranscriptomeBan Singleend \
        --outSAMattrRGline ID:$prefix SM:$prefix LB:library PL:illumina \
        --outFileNamePrefix $prefix 
        $args
    # Index the BAM file
    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
}