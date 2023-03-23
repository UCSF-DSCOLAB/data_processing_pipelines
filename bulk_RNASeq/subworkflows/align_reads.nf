//
// Read Alignment to reference genome
//

include { STAR_ALIGN } from '../modules/star_align'
include { PICARD_MARKDUPLICATES } from '../modules/picard_markduplicates'
include { SAMTOOLS_INDEX } from '../modules/samtools_index'
include { SAMTOOLS_SORT } from '../modules/samtools_sort'

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
    // Index BAM file
    //
    SAMTOOLS_INDEX (
        STAR_ALIGN.out.bam
    )
    //
    // Sort BAM file
    //
    SAMTOOLS_SORT (
        STAR_ALIGN.out.bam
    )
    

    emit:
    bam                 = SAMTOOLS_SORT.out.bam   
    bai                 = SAMTOOLS_INDEX.out.bai   
    transcriptome_bam   = STAR_ALIGN.out.transcriptome_bam
}