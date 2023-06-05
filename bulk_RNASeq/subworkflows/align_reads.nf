//
// Read Alignment to reference genome
//

include { STAR_ALIGN } from '../modules/star_align'
include { EXTRACT_MAPPED_READS as EXTRACT_MAPPED_GENOME_READS } from '../modules/samtools_extract_mapped_reads'
include { EXTRACT_MAPPED_READS as EXTRACT_MAPPED_TRANSCRIPTOME_READS } from '../modules/samtools_extract_mapped_reads'
include { SAMTOOLS_BAM_TO_CRAM as SAMTOOLS_GENOME_BAM_TO_CRAM } from '../modules/samtools_bam_to_cram'
include { SAMTOOLS_BAM_TO_CRAM as SAMTOOLS_TRANSCRIPTOME_BAM_TO_CRAM } from '../modules/samtools_bam_to_cram'
include { PICARD_MARKDUPLICATES } from '../modules/picard_markduplicates'
include { SAMTOOLS_INDEX } from '../modules/samtools_index'
include { SAMTOOLS_SORT } from '../modules/samtools_sort'
include { SAMTOOLS_STATS } from '../modules/samtools_stats'
include { SAMTOOLS_FLAGSTAT } from '../modules/samtools_flagstat'
include { SAMTOOLS_IDXSTATS } from '../modules/samtools_idxstats'

workflow ALIGN_READS {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    gtf                     // channel: /path/to/genome.gtf
    genome_dir

    main:
    
    //
    // Map reads with STAR
    //
    STAR_ALIGN (
        reads,
        gtf,
        genome_dir
    )
    //
    // Extract Mapped genome reads
    //
    EXTRACT_MAPPED_GENOME_READS (
        STAR_ALIGN.out.bam,
        prefix_addon = ''
    )
    //
    // Sort BAM file
    //
    SAMTOOLS_SORT (
        EXTRACT_MAPPED_GENOME_READS.out.mapped_bam
    )
    //
    // Index BAM file
    //
    SAMTOOLS_INDEX (
        EXTRACT_MAPPED_GENOME_READS.out.mapped_bam
    )
    //
    // Convert mapped genome reads from BAM to CRAM
    //
    SAMTOOLS_GENOME_BAM_TO_CRAM (
        EXTRACT_MAPPED_GENOME_READS.out.mapped_bam,
        params.genome,
        prefix_addon = ''
    )
    //
    // Extract Mapped transcriptome reads
    //
    EXTRACT_MAPPED_TRANSCRIPTOME_READS (
        STAR_ALIGN.out.transcriptome_bam,
        prefix_addon = '.transcriptome'
    )
    //
    // Convert mapped transcriptome reads from BAM to CRAM
    //
    SAMTOOLS_TRANSCRIPTOME_BAM_TO_CRAM (
        EXTRACT_MAPPED_TRANSCRIPTOME_READS.out.mapped_bam,
        params.transcript_fasta,
        prefix_addon = '.transcriptome'
    )
    //
    // Run samtools stats, flagstat and idxstats
    //
    ch_bam_bai = SAMTOOLS_SORT.out.bam
            .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
    SAMTOOLS_STATS ( ch_bam_bai, [] )
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    

    emit:
    bam                 = SAMTOOLS_SORT.out.bam   
    bai                 = SAMTOOLS_INDEX.out.bai   
    transcriptome_bam   = EXTRACT_MAPPED_TRANSCRIPTOME_READS.out.mapped_bam
    stats               = SAMTOOLS_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat            = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats            = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    log_final           = STAR_ALIGN.out.log_final

}