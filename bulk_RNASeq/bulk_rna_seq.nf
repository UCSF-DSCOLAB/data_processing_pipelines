#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Define the default parameters 
 */
params.input                    = ""
params.genome                   = ""
params.genome_idx               = ""
params.genome_dict              = ""
params.genome_dir               = ""
params.gtf                      = ""
params.transcript_fasta         = ""
params.gtf_group_features       = ""
params.gtf_extra_attributes     = ""
params.gatk_vf_cluster_size     = ""
params.gatk_vf_window_size      = ""
params.gatk_vf_fs_filter        = ""
params.gatk_vf_qd_filter        = ""
params.umitools_dedup_stats     = ""
params.dbsnp                    = ""
params.dbsnp_tbi                = ""
params.tmp_dir                  = ""
params.results_directory        = ""
params.rrna_db_file             = ""


// Check mandatory parameters (sample sheet)
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' } 

// Import SUBWORKFLOWS
include { INPUT_CHECK               } from './subworkflows/validate_input' 
include { ALIGN_READS               } from './subworkflows/align_reads'
include { BAM_DEDUP_UMITOOLS as BAM_DEDUP_UMITOOLS_GENOME } from './subworkflows/bam_dedup_umitools'
include { BAM_DEDUP_UMITOOLS as BAM_DEDUP_UMITOOLS_TRANSCRIPTOME } from './subworkflows/bam_dedup_umitools'
include { BAM_MARKDUPLICATES_PICARD } from './subworkflows/post_process_bam'
include { QUANTIFY_SALMON           } from './subworkflows/quantify_transcriptome'

// Import MODULES
include { CAT_FASTQ                 } from './modules/cat_fastq'
include { FASTP_TRIM_ADAPTERS       } from './modules/fastp_trim_adapters'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANSCRIPTOME_PRE_DEDUP } from './modules/samtools_sort'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANSCRIPTOME_POST_DEDUP } from './modules/samtools_sort'
include { SAMTOOLS_INDEX            } from './modules/samtools_index'
include { UMITOOLS_PREPARE_FOR_SALMON } from './modules/umitools_prepare_for_salmon'
include { SORTMERNA                 } from './modules/sortmerna_rrna_removal'
include { KALLISTO_QUANT            } from './modules/kallisto_quant'
include { GATK4_SPLITNCIGARREADS    } from './modules/gatk4_splitncigar'
include { GATK4_BASE_RECALIBRATOR   } from './modules/gatk4_recalibrator'
include { GATK4_APPLY_BQSR          } from './modules/gatk4_apply_bqsr'
include { GATK4_HAPLOTYPECALLER     } from './modules/gatk4_haplotype_caller'
include { GATK4_VARIANTFILTRATION   } from './modules/gatk4_variant_filter'
include { MULTIQC                   } from './modules/multiqc'


workflow {
    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    //
    // MODULE: Trim adapter sequences from FastQ reads
    //
    ch_trimmed_reads = Channel.empty()
    ch_trim_multiqc = Channel.empty()
    FASTP_TRIM_ADAPTERS (
        ch_cat_fastq
    )
    ch_trimmed_reads = FASTP_TRIM_ADAPTERS.out.trimmed_reads
    ch_trim_multiqc = FASTP_TRIM_ADAPTERS.out.log
    ch_reports = ch_reports.mix(FASTP_TRIM_ADAPTERS.out.log.collect{it[1]}.ifEmpty([]))
    //
    // MODULE: Remove ribosomal RNA reads
    //
    // ch_filtered_reads = Channel.empty()
    ch_sortmerna_multiqc = Channel.empty()
    if (params.rrna_db_file && params.filter_rrna) {
        sortmerna_fastas_data = file(params.rrna_db_file).readLines()
        lst_sortmerna_fastas = sortmerna_fastas_data.collect { file(it) }
        SORTMERNA (
            ch_trimmed_reads,
            lst_sortmerna_fastas
        )
        ch_trimmed_reads = SORTMERNA.out.reads
        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_reports = ch_reports.mix(SORTMERNA.out.log_final.collect{it[1]}.ifEmpty([]))
    }
    //
    // MODULE: Quantify transcriptome abundance using Kallisto
    //
    KALLISTO_QUANT(
        ch_trimmed_reads
    )
    // SUBWORKFLOW: Align FastQ reads; sort, and index BAM files
    //
    ch_star_bam = Channel.empty()
    ch_star_bai = Channel.empty()
    ch_star_stats    = Channel.empty()
    ch_star_flagstat = Channel.empty()
    ch_star_idxstats = Channel.empty()
    ch_star_multiqc  = Channel.empty()
    ch_transcriptome_bam_star = Channel.empty()
    ALIGN_READS(
        ch_trimmed_reads,
        params.gtf,
        params.genome_dir
    )
    ch_star_bam = ALIGN_READS.out.bam
    ch_star_bai = ALIGN_READS.out.bai
    ch_transcriptome_bam_star = ALIGN_READS.out.transcriptome_bam
    ch_star_stats    = ALIGN_READS.out.stats
    ch_star_flagstat = ALIGN_READS.out.flagstat
    ch_star_idxstats = ALIGN_READS.out.idxstats
    ch_star_multiqc  = ALIGN_STAR.out.log_final
    ch_reports = ch_reports.mix(ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
    ch_star_bam_bai = ch_star_bam.join(ch_star_bai, by: [0])
    // Deduplicate genome BAM file before downstream analysis
    // BAM_DEDUP_UMITOOLS_GENOME (
    //     ch_star_bam_bai,
    //     params.umitools_dedup_stats
    // )
    // ch_genome_bam        = BAM_DEDUP_UMITOOLS_GENOME.out.bam
    // ch_genome_bam_index  = BAM_DEDUP_UMITOOLS_GENOME.out.bai
    // ch_samtools_stats    = BAM_DEDUP_UMITOOLS_GENOME.out.stats
    // ch_samtools_flagstat = BAM_DEDUP_UMITOOLS_GENOME.out.flagstat
    // ch_samtools_idxstats = BAM_DEDUP_UMITOOLS_GENOME.out.idxstats
    // Co-ordinate sort, index and run stats on transcriptome BAM
    ch_transcriptome_bam_sort = Channel.empty()
    SAMTOOLS_SORT_TRANSCRIPTOME_PRE_DEDUP ( ch_transcriptome_bam_star )
    ch_transcriptome_bam_sort = SAMTOOLS_SORT_TRANSCRIPTOME_PRE_DEDUP.out.bam
    SAMTOOLS_INDEX ( ch_transcriptome_bam_sort )
    // ch_sam_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    ch_transcriptome_bam_sort
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
                }
        .set { ch_transcriptome_bam_pre }
    // ch_transcriptome_sorted_bam = SAMTOOLS_SORT_TRANSCRIPTOME_PRE_DEDUP.out.bam
    // ch_transcriptome_sorted_bai = SAMTOOLS_INDEX.out.bai
    // ch_transcriptome_sorted_bam_bai = ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0])
    // // Deduplicate transcriptome BAM file before read counting with Salmon
    // BAM_DEDUP_UMITOOLS_TRANSCRIPTOME (
    //     ch_transcriptome_sorted_bam_bai,
    //     params.umitools_dedup_stats
    // )
    // // Name sort BAM before passing to Salmon
    // SAMTOOLS_SORT_TRANSCRIPTOME_POST_DEDUP (
    //     ch_transcriptome_sorted_bam
    // )
    // Only run prepare_for_rsem.py on paired-end BAM files
    // SAMTOOLS_SORT_TRANSCRIPTOME_POST_DEDUP
    //     .out
    //     .bam
    //     .branch {
    //         meta, bam ->
    //             single_end: meta.single_end
    //                 return [ meta, bam ]
    //             paired_end: !meta.single_end
    //                 return [ meta, bam ]
    //     }
    //     .set { ch_umitools_dedup_bam }
    // Fix paired-end reads in name sorted BAM file
    // See: https://github.com/nf-core/rnaseq/issues/828
    UMITOOLS_PREPARE_FOR_SALMON (
        ch_transcriptome_bam_pre.paired_end
    )
    ch_transcriptome_bam_pre
        .single_end
        .mix(UMITOOLS_PREPARE_FOR_SALMON.out.bam)
        .set { ch_transcriptome_bam }
    //
    // SUBWORKFLOW: Count reads from BAM alignments using Salmon
    //
    // QUANTIFY_SALMON (
    //     ch_transcriptome_bam,
    //     params.transcript_fasta,
    //     params.gtf,
    //     params.salmon_quant_libtype ?: ''
    // )
    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    BAM_MARKDUPLICATES_PICARD (
        ch_star_bam_bai,
        params.genome,
        params.genome_idx
    )
    ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
    ch_genome_bai             = BAM_MARKDUPLICATES_PICARD.out.bai
    ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
    ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
    ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
    ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
    ch_reports = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]))
    ch_reports = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_genome_bam_bai = ch_genome_bam.join(ch_genome_bai, by: [0])
    //
    // MODULE: SplitNCigarReads and reassign mapping qualities
    //
    ch_split_bam = Channel.empty()
    ch_split_bai = Channel.empty()
    GATK4_SPLITNCIGARREADS (
        ch_genome_bam_bai,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_split_bam = GATK4_SPLITNCIGARREADS.out.bam
    ch_split_bai = GATK4_SPLITNCIGARREADS.out.bai
    //
    // MODULE: Base Recalibration table generation
    //
    ch_recal_table = Channel.empty()
    GATK4_BASE_RECALIBRATOR (
        ch_split_bam,
        ch_split_bai,
        params.genome,
        params.genome_idx,
        params.genome_dict,
        params.dbsnp,
        params.dbsnp_tbi
    )
    ch_recal_table = GATK4_BASE_RECALIBRATOR.out.table
    ch_reports = ch_reports.mix(ch_recal_table.map{ meta, table -> table})
    //
    // MODULE: Apply BQSR using recalibration table
    //
    ch_split_bam_bai = ch_split_bam.join(ch_split_bai, by: [0])
    ch_bam_bai_bqsr = ch_split_bam_bai.join(ch_recal_table, by: [0])
    ch_bam_variant_calling = Channel.empty()
    ch_bai_variant_calling = Channel.empty()
    GATK4_APPLY_BQSR (
        ch_bam_bai_bqsr,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_bam_variant_calling = GATK4_APPLY_BQSR.out.bam
    ch_bai_variant_calling = GATK4_APPLY_BQSR.out.bai
    ch_reports = ch_reports.mix(GATK4_APPLY_BQSR.out.qc.collect{it[1]}.ifEmpty([]))
    //
    // MODULE: Call SNPs and Indels using HaplotypeCaller
    //
    ch_bam_bai_variant_calling = ch_bam_variant_calling.join(ch_bai_variant_calling, by: [0])
    ch_haplotype_vcf = Channel.empty()
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_variant_calling,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_haplotype_vcf = GATK4_HAPLOTYPECALLER.out.vcf
    ch_haplotype_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    ch_haplotype_vcf_tbi = ch_haplotype_vcf.join(ch_haplotype_tbi, by: [0])
    //
    // MODULE: Filter variants using VariantFiltration
    //
    ch_filtered_vcf = Channel.empty()
    GATK4_VARIANTFILTRATION (
        ch_haplotype_vcf_tbi,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf 
    //
    // MODULE: Generate QC reports using MULTIQC
    //
    ch_multiqc_files = Channel
                            .empty()
                            .mix(ch_reports.collect())
    MULTIQC (ch_multiqc_files.collect())
    multiqc_report = MULTIQC.out.report.toList()
}