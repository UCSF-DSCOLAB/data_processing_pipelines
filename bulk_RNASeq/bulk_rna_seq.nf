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
include { BAM_MARKDUPLICATES_PICARD } from './subworkflows/post_process_bam'
include { QUANTIFY_SALMON           } from './subworkflows/quantify_transcriptome'

// Import MODULES
include { CAT_FASTQ                 } from './modules/cat_fastq'
include { FASTP_TRIM_ADAPTERS       } from './modules/fastp_trim_adapters'
include { SORTMERNA                 } from './modules/sortmerna_rrna_removal'
include { GATK4_SPLITNCIGARREADS    } from './modules/gatk4_splitncigar'
include { GATK4_BASE_RECALIBRATOR   } from './modules/gatk4_recalibrator'
include { GATK4_APPLY_BQSR          } from './modules/gatk4_apply_bqsr'
include { GATK4_HAPLOTYPECALLER     } from './modules/gatk4_haplotype_caller'
include { GATK4_VARIANTFILTRATION   } from './modules/gatk4_variant_filter'


workflow {
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
    FASTP_TRIM_ADAPTERS (
        ch_cat_fastq
    )
    ch_trimmed_reads = FASTP_TRIM_ADAPTERS.out.trimmed_reads
    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()
    if (params.rrna_db_file) {
        sortmerna_fastas_data = file(params.rrna_db_file).readLines()
        lst_sortmerna_fastas = sortmerna_fastas_data.collect { file(it) }

        SORTMERNA (
            ch_trimmed_reads,
            lst_sortmerna_fastas
        )
        .reads
        .set { ch_filtered_reads }

        ch_sortmerna_multiqc = SORTMERNA.out.log
    }
    //
    // SUBWORKFLOW: Align FastQ reads; sort, and index BAM files
    //
    ch_genome_bam = Channel.empty()
    ch_genome_bai = Channel.empty()
    ALIGN_READS(
        ch_filtered_reads,
        params.gtf,
        params.genome_dir
    )
    ch_star_bam = ALIGN_READS.out.bam
    ch_star_bai = ALIGN_READS.out.bai
    ch_transcriptome_bam = ALIGN_READS.out.transcriptome_bam
    ch_star_bam_bai = ch_star_bam.join(ch_star_bai, by: [0])
    //
    // SUBWORKFLOW: Count reads from BAM alignments using Salmon
    //
    QUANTIFY_SALMON (
        ch_transcriptome_bam,
        params.transcript_fasta,
        params.gtf,
        params.salmon_quant_libtype ?: ''
    )
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
}