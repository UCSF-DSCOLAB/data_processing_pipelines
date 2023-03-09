//
// Read Alignment to reference genome
//

include { STAR_ALIGN } from '../modules/star_align'

workflow ALIGN_READS {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    gtf                     // channel: /path/to/genome.gtf
    genome_dir
    tmp_dir

    main:

    //
    // Map reads with STAR
    //
    STAR_ALIGN (
        reads,
        gtf,
        genome_dir,
        tmp_dir
    )

    emit:
    bam              = STAR_ALIGN.out.bam             // channel: [ val(meta), bam            ]
    bai              = STAR_ALIGN.out.bai             // channel: [ val(meta), bai            ]

}