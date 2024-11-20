process STAR_ALIGN {
    tag "$meta.id"
    clusterOptions = '-S /bin/bash'
    label 'star_align'
    memory {
        if (meta.single_end) {
          // File size in GB
          fileSize = reads.size() / (1024 * 1024 * 1024)
        } else {
          // File size in GB
          fileSize = reads[0].size() / (1024 * 1024 * 1024)
        }
	return 20.GB * (1 + (fileSize*0.1))
    }
    publishDir "${params.results_directory}/star", mode: 'copy', pattern: "${prefix}ReadsPerGene.out.tab"
    publishDir "${params.results_directory}/star", mode: 'copy', pattern: "${prefix}Log.final.out"

    input:
    tuple val(meta), path(reads)
    path gtf
    path genomeDir

    output:
    tuple val(meta), path('*sortedByCoord.out.bam'), emit: bam
    tuple val(meta), path('*toTranscriptome.out.bam'), emit: transcriptome_bam
    tuple val(meta), path('*ReadsPerGene.out.tab'), emit: gene_counts 
    tuple val(meta), path('*Log.final.out'), emit: log_final


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \\
        --readFilesIn $reads  \\
        --genomeDir $genomeDir \\
        --runThreadN $task.cpus \\
        --sjdbGTFfile $gtf \\
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outReadsUnmapped None \
	    --outSAMunmapped Within KeepPairs \
        --outSAMattrRGline ID:$prefix SM:$prefix LB:library PL:illumina \
        --outFileNamePrefix $prefix 
        $args
    """
}