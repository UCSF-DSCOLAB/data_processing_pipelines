#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Display the default parameters (configure via nextflow.config)
 */
params.input                    = ""
params.genome                   = ""
params.genome_idx               = ""
params.genome_dict              = ""
params.genome_dir               = ""
params.gtf                      = ""
params.transcript_fasta         = ""
params.transcript_index         = ""
params.gtf_group_features       = ""
params.gtf_extra_attributes     = ""
params.fragment_length_mean     = ""
params.fragment_length_std      = ""
params.gatk_vf_cluster_size     = ""
params.gatk_vf_window_size      = ""
params.gatk_vf_fs_filter        = ""
params.gatk_vf_qd_filter        = ""
params.umitools_dedup_stats     = ""
params.dbsnp                    = ""
params.dbsnp_tbi                = ""
params.contig_format_map        = ""
params.format_contigs           = ""
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
include { KALLISTO_QUANT            } from './modules/kallisto_quant'
include { GATK4_SPLITNCIGARREADS    } from './modules/gatk4_splitncigar'
include { GATK4_BASE_RECALIBRATOR   } from './modules/gatk4_recalibrator'
include { GATK4_APPLY_BQSR          } from './modules/gatk4_apply_bqsr'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR } from './modules/samtools_index'
include { GATK4_HAPLOTYPECALLER     } from './modules/gatk4_haplotype_caller'
include { GATK4_VARIANTFILTRATION   } from './modules/gatk4_variant_filter'
include { BCFTOOLS_CONTIG_CONVERSION} from './modules/bcftools_contig_conversion'
include { BCFTOOLS_SORT_VCF   }       from './modules/bcftools_sort_vcf'
include { BCFTOOLS_INDEX_VCF   }      from './modules/bcftools_index_vcf'
include { BCFTOOLS_MERGE_VCF        } from './modules/bcftools_merge_vcf'
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
    ch_trim_multiqc = FASTP_TRIM_ADAPTERS.out.json_report
    ch_reports = ch_reports.mix(ch_trim_multiqc)
    //
    // MODULE: Remove ribosomal RNA reads
    //
    if (params.rrna_db_file && params.filter_rrna) {
        ch_sortmerna_multiqc = Channel.empty()
        sortmerna_fastas_data = file(params.rrna_db_file).readLines()
        lst_sortmerna_fastas = sortmerna_fastas_data.collect { file(it) }
        SORTMERNA (
            ch_trimmed_reads,
            lst_sortmerna_fastas
        )
        ch_trimmed_reads = SORTMERNA.out.reads
        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_reports = ch_reports.mix(SORTMERNA.out.log.collect{it[1]}.ifEmpty([]))
    }
    //
    // MODULE: Quantify transcriptome abundance using Kallisto
    //
    ch_kallisto_multiqc = Channel.empty()
    KALLISTO_QUANT(
        ch_trimmed_reads
    )
    ch_kallisto_multiqc = KALLISTO_QUANT.out.log
    ch_reports = ch_reports.mix(KALLISTO_QUANT.out.log.collect{it[1]}.ifEmpty([]))
    //
    // SUBWORKFLOW: Align FastQ reads; sort, and index BAM files
    //
    ch_star_bam = Channel.empty()
    ch_star_bai = Channel.empty()
    ch_transcriptome_bam_star = Channel.empty()
    ch_star_stats    = Channel.empty()
    ch_star_flagstat = Channel.empty()
    ch_star_idxstats = Channel.empty()
    ch_star_multiqc  = Channel.empty()
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
    ch_star_multiqc  = ALIGN_READS.out.log_final
    ch_reports = ch_reports.mix(ALIGN_READS.out.log_final.collect{it[1]}.ifEmpty([]))
    ch_star_bam_bai = ch_star_bam.join(ch_star_bai, by: [0])
    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_genome_bam             = Channel.empty()
    ch_genome_bai             = Channel.empty()
    ch_samtools_stats         = Channel.empty()
    ch_samtools_flagstat      = Channel.empty()
    ch_samtools_idxstats      = Channel.empty()
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
    // MODULE: Apply BQSR using recalibration table, then index
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
    SAMTOOLS_INDEX_BQSR (
        GATK4_APPLY_BQSR.out.bam
    )
    ch_bam_variant_calling = GATK4_APPLY_BQSR.out.bam
    ch_bai_variant_calling = SAMTOOLS_INDEX_BQSR.out.bai
    // ch_reports = ch_reports.mix(GATK4_APPLY_BQSR.out.qc.collect{it[1]}.ifEmpty([]))
    //
    // MODULE: Call SNPs and Indels using HaplotypeCaller
    //
    ch_bam_bai_variant_calling = ch_bam_variant_calling.join(ch_bai_variant_calling, by: [0])
    ch_haplotype_vcf = Channel.empty()
    ch_haplotype_tbi = Channel.empty()
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_variant_calling,
        params.genome,
        params.genome_idx,
        params.genome_dict,
        params.dbsnp,
        params.dbsnp_tbi
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
    if (params.format_contigs && params.contig_format_map) {
        //
        // MODULE: Convert VCF contigs to desired naming format (e.g. ucsc)
        //
        BCFTOOLS_CONTIG_CONVERSION (
           ch_filtered_vcf
        )
        ch_filtered_vcf = BCFTOOLS_CONTIG_CONVERSION.out.formatted_vcf
    }
    //
    // MODULE: Sort and index VCFs
    //
    ch_sorted_filtered_vcf = Channel.empty()
    BCFTOOLS_SORT_VCF (
        ch_filtered_vcf
    )
    ch_sorted_vcf = BCFTOOLS_SORT_VCF.out.sorted_vcf
    //
    // MODULE: Index VCFs
    //
    ch_vcf_index = Channel.empty()
    BCFTOOLS_INDEX_VCF (
        ch_sorted_vcf
    )
    ch_vcf_index = BCFTOOLS_INDEX_VCF.out.vcf_index
    ch_vcf = ch_sorted_vcf.join(ch_vcf_index, by: [0])
    // Collect all VCFs and index files from upstream process
    meta = ch_vcf
    .map { tuple -> tuple[0]}
    .collect()
    vcfs = ch_vcf
    .map { tuple -> tuple[1]}
    .collect()
    tbis = ch_vcf
    .map { tuple -> tuple[2]}
    .collect()
    //
    // MODULE: Merge VCFs
    //
    BCFTOOLS_MERGE_VCF (
        meta, 
        vcfs, 
        tbis
    )
    //
    // MODULE: Generate QC reports using MULTIQC
    //
    ch_multiqc_files = Channel
                            .empty()
                            .mix(ch_reports.collect())
    MULTIQC (ch_multiqc_files.collect())
    multiqc_report = MULTIQC.out.report.toList()
}