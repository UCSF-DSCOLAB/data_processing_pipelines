process STAR_ALIGN {
    tag "$meta.id"
    // clusterOptions = '-S /bin/bash'
    label 'star_align', 'per_sample'
    memory {
      /* Size of the reads in GiB */
      def fileSizeGB = meta.single_end
                      ? (reads.size() / (1024 ** 3))
                      : (reads[0].size()    / (1024 ** 3))

      /* Formula-based requirement: 25 GB × ( 2 + 0.1 × size ) */
      def required = 25.GB * ( 2 + (fileSizeGB * 0.1) )

      /* Never ask for less than 8 GB */
      return [required, 25.GB].max()          // equivalent to Math.max(required, 8.GB)
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
        --outFileNamePrefix $prefix \\
        --outFilterMismatchNoverLmax ${params.star_outfilter_mismatch_n_over_lmax} \\
        --alignSJoverhangMin ${params.star_align_sjoverhang_min} \\
        --outFilterMultimapNmax ${params.star_outfilter_multimap_nmax} \\
        --seedSearchStartLmax ${params.star_seed_search_start_lmax} \\
        ${params.star_additional ?: ''} \\
        $args
    """
}