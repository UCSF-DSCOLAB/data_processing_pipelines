//
// Pseudo-alignment and quantification with Salmon
//

include { SALMON_QUANT    } from '../modules/salmon_quant'
include { SALMON_TX2GENE  } from '../modules/salmon_tx2gene'
include { SALMON_TXIMPORT } from '../modules/salmon_tximport'

include { SALMON_SUMMARIZE_EXPERIMENT as SALMON_SE_GENE               } from '../modules/salmon_summarize'
include { SALMON_SUMMARIZE_EXPERIMENT as SALMON_SE_GENE_LENGTH_SCALED } from '../modules/salmon_summarize'
include { SALMON_SUMMARIZE_EXPERIMENT as SALMON_SE_GENE_SCALED        } from '../modules/salmon_summarize'
include { SALMON_SUMMARIZE_EXPERIMENT as SALMON_SE_TRANSCRIPT         } from '../modules/salmon_summarize'

workflow QUANTIFY_SALMON {
    take:
    bam            // channel: [ val(meta), [ bam ] ]
    transcript_fasta // channel: /path/to/transcript.fasta
    gtf              // channel: /path/to/genome.gtf
    lib_type         //     val: String to override salmon library type

    main:
    //
    // Quantify and merge counts across samples
    //
    SALMON_QUANT ( bam, 
                   gtf, 
                   transcript_fasta,  
                   lib_type 
    )

    SALMON_TX2GENE ( SALMON_QUANT.out.results.collect{it[1]}, gtf )

    SALMON_TXIMPORT ( SALMON_QUANT.out.results.collect{it[1]}, SALMON_TX2GENE.out.tsv.collect() )

    SALMON_SE_GENE (
        SALMON_TXIMPORT.out.counts_gene,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )

    SALMON_SE_GENE_LENGTH_SCALED (
        SALMON_TXIMPORT.out.counts_gene_length_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )

    SALMON_SE_GENE_SCALED (
        SALMON_TXIMPORT.out.counts_gene_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )

    SALMON_SE_TRANSCRIPT (
        SALMON_TXIMPORT.out.counts_transcript,
        SALMON_TXIMPORT.out.tpm_transcript,
        SALMON_TX2GENE.out.tsv.collect()
    )

    emit:
    results                       = SALMON_QUANT.out.results                      // channel: [ val(meta), results_dir ]

    tpm_gene                      = SALMON_TXIMPORT.out.tpm_gene                  // channel: [ val(meta), counts ]
    counts_gene                   = SALMON_TXIMPORT.out.counts_gene               // channel: [ val(meta), counts ]
    counts_gene_length_scaled     = SALMON_TXIMPORT.out.counts_gene_length_scaled // channel: [ val(meta), counts ]
    counts_gene_scaled            = SALMON_TXIMPORT.out.counts_gene_scaled        // channel: [ val(meta), counts ]
    tpm_transcript                = SALMON_TXIMPORT.out.tpm_transcript            // channel: [ val(meta), counts ]
    counts_transcript             = SALMON_TXIMPORT.out.counts_transcript         // channel: [ val(meta), counts ]

    merged_gene_rds               = SALMON_SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled = SALMON_SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled        = SALMON_SE_GENE_SCALED.out.rds                 //    path: *.rds

    merged_counts_transcript      = SALMON_TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = SALMON_TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds         = SALMON_SE_TRANSCRIPT.out.rds                  //    path: *.rds
}
