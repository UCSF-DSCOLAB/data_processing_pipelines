//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../modules/picard_markduplicates'
include { SAMTOOLS_INDEX        } from '../modules/samtools_index'
include { SAMTOOLS_STATS        } from '../modules/samtools_stats'
include { SAMTOOLS_FLAGSTAT     } from '../modules/samtools_flagstat'
include { SAMTOOLS_IDXSTATS     } from '../modules/samtools_idxstats'

workflow BAM_MARKDUPLICATES_PICARD {

    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ fasta ]
    ch_fai   // channel: [ fai ]

    main:
    //
    // MODULE: Mark duplicate reads using Picard
    //
    PICARD_MARKDUPLICATES ( ch_bam, ch_fasta, ch_fai )
    //
    // MODULE: Index reads using samtools
    //
    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    PICARD_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }
    //
    // MODULE: Generate stats, flag, and index them using samtools
    //
    SAMTOOLS_STATS ( ch_bam_bai, ch_fasta )
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    SAMTOOLS_IDXSTATS ( ch_bam_bai )

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), [ bam ] ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]

    stats    = SAMTOOLS_STATS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
}