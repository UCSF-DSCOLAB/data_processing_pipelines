#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 * Display the default parameters (configure via nextflow.config)
*/
params.input_sample_sheet       = ""
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
params.gene_mapper              = ""
params.contig_format_map        = ""
params.format_contigs           = ""
params.tmp_dir                  = ""
params.results_directory        = ""
params.rrna_db_file             = ""

// Check mandatory parameters (sample sheet)
if (params.input_sample_sheet) { ch_input = file(params.input_sample_sheet) } else { exit 1, 'Input samplesheet not specified!' } 
if ((params.demux_snp_panel && !params.demux_snp_panel_tbi) ||
    (!params.demux_snp_panel && params.demux_snp_panel_tbi)) {
    exit 1, 'demux_snp_panel and demux_snp_panel_tbi must be supplied together!'
}
if (params.demux_max_missing_fraction < 0 || params.demux_max_missing_fraction > 1) {
    exit 1, 'demux_max_missing_fraction must be between 0 and 1!'
}
//if (params.transcript_index) { ch_transcript_index = file(params.transcript_index) } else { exit 1, 'Input transcript index not specified!' } 

// Import SUBWORKFLOWS
include { INPUT_CHECK               } from './subworkflows/validate_input' 
include { ALIGN_READS               } from './subworkflows/align_reads'
include { BAM_MARKDUPLICATES_PICARD } from './subworkflows/post_process_bam'
// include { QUANTIFY_SALMON           } from './subworkflows/quantify_transcriptome'

// Import MODULES
include { CAT_FASTQ                 } from './modules/cat_fastq'
include { FASTP_TRIM_ADAPTERS       } from './modules/fastp_trim_adapters'
include { FASTP_PREPARE_GENOTYPING  } from './modules/fastp_prepare_genotyping'
include { SORTMERNA_RIBOSOMAL_RNA_REMOVAL       } from './modules/sortmerna_rrna_removal'
include { KALLISTO_QUANT            } from './modules/kallisto_quant'
include { CUSTOM_MERGE_COUNTS       } from './modules/custom_merge_counts'
include { GATK4_SPLITNCIGARREADS    } from './modules/gatk4_splitncigar'
include { GATK4_BASE_RECALIBRATOR   } from './modules/gatk4_recalibrator'
include { GATK4_APPLY_BQSR          } from './modules/gatk4_apply_bqsr'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR } from './modules/samtools_index'
include { GATK4_HAPLOTYPECALLER     } from './modules/gatk4_haplotype_caller'
include { GATK4_COMBINEGVCFS        } from './modules/gatk4_combine_gvcfs'
include { GATK4_GENOTYPEGVCFS       } from './modules/gatk4_genotype_gvcfs'
include { GATK4_SELECT_SNPS         } from './modules/gatk4_select_snps'
include { GATK4_VARIANTFILTRATION   } from './modules/gatk4_variant_filter'
include { GATK4_SELECT_PASS_VARIANTS} from './modules/gatk4_select_pass_variants'
include { BCFTOOLS_DEMUX_FILTER     } from './modules/bcftools_demux_filter'
include { BCFTOOLS_INTERSECT_PANEL  } from './modules/bcftools_intersect_panel'
include { BCFTOOLS_EXCLUDE_REGIONS  } from './modules/bcftools_exclude_regions'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_RAW   } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FINAL } from './modules/bcftools_stats'
include { GENOTYPE_VCF_QC           } from './modules/genotype_vcf_qc'
include { BCFTOOLS_CONTIG_CONVERSION} from './modules/bcftools_contig_conversion'
include { BCFTOOLS_SORT_VCF   }       from './modules/bcftools_sort_vcf'
include { BCFTOOLS_INDEX_VCF   }      from './modules/bcftools_index_vcf'
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
            meta_clone.id = meta_clone.id.replaceFirst(/_T\d+$/, '')
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
    // MODULE: Prepare an independent genotype-safe read branch. Overlap
    // correction is intentionally omitted so observed alleles are not rewritten.
    //
    FASTP_PREPARE_GENOTYPING (ch_cat_fastq)
    ch_genotyping_reads = FASTP_PREPARE_GENOTYPING.out.trimmed_reads
    ch_reports = ch_reports.mix(FASTP_PREPARE_GENOTYPING.out.json_report)
    //
    // MODULE: Remove ribosomal RNA reads
    //
    if (params.rrna_db_file && params.filter_rrna) {
        ch_sortmerna_multiqc = Channel.empty()
        sortmerna_fastas_data = file(params.rrna_db_file).readLines()
        lst_sortmerna_fastas = sortmerna_fastas_data.collect { file(it) }
        SORTMERNA_RIBOSOMAL_RNA_REMOVAL (
            ch_trimmed_reads,
            lst_sortmerna_fastas
        )
        ch_trimmed_reads = SORTMERNA_RIBOSOMAL_RNA_REMOVAL.out.reads
        ch_sortmerna_multiqc = SORTMERNA_RIBOSOMAL_RNA_REMOVAL.out.log
        ch_reports = ch_reports.mix(SORTMERNA_RIBOSOMAL_RNA_REMOVAL.out.log.map{it[1]}.ifEmpty([]))
    }
    //
    // MODULE: Quantify transcriptome abundance using Kallisto
    //
    ch_kallisto_multiqc = Channel.empty()
    ch_kallisto_counts = Channel.empty()
    KALLISTO_QUANT(
        ch_trimmed_reads,
        params.transcript_index
    )
    ch_kallisto_counts = KALLISTO_QUANT.out.abundance_tsv
    ch_kallisto_multiqc = KALLISTO_QUANT.out.log
    ch_reports = ch_reports.mix(KALLISTO_QUANT.out.log.map{it[1]}.ifEmpty([]))
    //
    // MODULE: Merge all transcriptome quantification into a single file
    //
    // counts = ch_kallisto_counts
    //     .filter { it[1]?.exists() }
    //     .map { tuple -> tuple[1] }

    CUSTOM_MERGE_COUNTS(
        ch_kallisto_counts.map { tuple -> tuple[1] }
    )
    // counts = ch_kallisto_counts.map { tuple -> tuple[1]}.collect()
    // CUSTOM_MERGE_COUNTS (
    //     counts
    // )
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
        ch_genotyping_reads,
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
    ch_reports = ch_reports.mix(ALIGN_READS.out.log_final.map{it[1]}.ifEmpty([]))
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
    ch_reports = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.stats.map{it[1]}.ifEmpty([]))
    ch_reports = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.map{it[1]}.ifEmpty([]))
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
    ch_bam_variant_calling = Channel.empty()
    ch_bai_variant_calling = Channel.empty()
    if (params.run_bqsr) {
        //
        // MODULES: Generate and apply BQSR. This remains configurable so its
        // effect on orthogonal genotype concordance can be benchmarked.
        //
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
        ch_split_bam_bai = ch_split_bam.join(ch_split_bai, by: [0])
        ch_bam_bai_bqsr = ch_split_bam_bai.join(ch_recal_table, by: [0])
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
    } else {
        ch_bam_variant_calling = ch_split_bam
        ch_bai_variant_calling = ch_split_bai
    }
    //
    // MODULE: Emit a reference-confidence GVCF for each donor.
    //
    ch_bam_bai_variant_calling = ch_bam_variant_calling.join(ch_bai_variant_calling, by: [0])
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_variant_calling,
        params.genome,
        params.genome_idx,
        params.genome_dict,
        params.dbsnp,
        params.dbsnp_tbi
    )
    ch_haplotype_gvcf = GATK4_HAPLOTYPECALLER.out.vcf
    ch_haplotype_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    ch_haplotype_gvcf_tbi = ch_haplotype_gvcf.join(ch_haplotype_tbi, by: [0])

    //
    // MODULES: Combine reference-confidence records and jointly genotype donors.
    // This replaces bcftools merge of independently discovered variant-only VCFs.
    //
    cohort_meta = Channel.value([id: 'cohort'])
    cohort_gvcfs = ch_haplotype_gvcf_tbi.map { meta, gvcf, tbi -> gvcf }.collect()
    cohort_gvcf_tbis = ch_haplotype_gvcf_tbi.map { meta, gvcf, tbi -> tbi }.collect()
    GATK4_COMBINEGVCFS (
        cohort_meta,
        cohort_gvcfs,
        cohort_gvcf_tbis,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    GATK4_GENOTYPEGVCFS (
        GATK4_COMBINEGVCFS.out.gvcf,
        params.genome,
        params.genome_idx,
        params.genome_dict,
        params.dbsnp,
        params.dbsnp_tbi
    )

    //
    // MODULES: retain biallelic SNPs, mask low-quality donor genotypes, and
    // physically select PASS sites.
    //
    GATK4_SELECT_SNPS (
        GATK4_GENOTYPEGVCFS.out.vcf,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    GATK4_VARIANTFILTRATION (
        GATK4_SELECT_SNPS.out.vcf,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    GATK4_SELECT_PASS_VARIANTS (
        GATK4_VARIANTFILTRATION.out.vcf,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )

    //
    // MODULES: normalize, remove cohort-monomorphic/high-missingness sites,
    // optionally intersect a common germline panel, and exclude difficult regions.
    //
    BCFTOOLS_DEMUX_FILTER (
        GATK4_SELECT_PASS_VARIANTS.out.vcf,
        params.genome,
        params.genome_idx
    )
    ch_demux_ready_vcf = BCFTOOLS_DEMUX_FILTER.out.vcf

    if (params.demux_snp_panel && params.demux_snp_panel_tbi) {
        BCFTOOLS_INTERSECT_PANEL (
            ch_demux_ready_vcf,
            params.demux_snp_panel,
            params.demux_snp_panel_tbi,
            params.genome,
            params.genome_idx
        )
        ch_demux_ready_vcf = BCFTOOLS_INTERSECT_PANEL.out.vcf
    }

    if (params.demux_exclude_regions) {
        BCFTOOLS_EXCLUDE_REGIONS (
            ch_demux_ready_vcf,
            params.demux_exclude_regions
        )
        ch_demux_ready_vcf = BCFTOOLS_EXCLUDE_REGIONS.out.vcf
    }

    // Rename contigs only after all reference-based normalization/filtering.
    if (params.format_contigs && params.contig_format_map) {
        BCFTOOLS_CONTIG_CONVERSION (
            ch_demux_ready_vcf.map { meta, vcf, tbi -> [meta, vcf] },
            params.contig_format_map
        )
        BCFTOOLS_SORT_VCF (BCFTOOLS_CONTIG_CONVERSION.out.formatted_vcf)
        BCFTOOLS_INDEX_VCF (BCFTOOLS_SORT_VCF.out.sorted_vcf)
    }

    //
    // MODULES: raw and final cohort/genotype QC for MultiQC and review.
    //
    BCFTOOLS_STATS_RAW (
        GATK4_GENOTYPEGVCFS.out.vcf.map { meta, vcf, tbi -> ['raw_joint', vcf, tbi] }
    )
    BCFTOOLS_STATS_FINAL (
        ch_demux_ready_vcf.map { meta, vcf, tbi -> ['demux_ready', vcf, tbi] }
    )
    GENOTYPE_VCF_QC (ch_demux_ready_vcf)
    ch_reports = ch_reports.mix(BCFTOOLS_STATS_RAW.out.stats.map { label, report -> report })
    ch_reports = ch_reports.mix(BCFTOOLS_STATS_FINAL.out.stats.map { label, report -> report })
    //
    // MODULE: Generate QC reports using MULTIQC
    //
    // After correcting all instances, you can now filter and use ch_reports
    // ch_multiqc_files = ch_reports.filter { it.exists() }
    // MULTIQC(ch_reports)
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_multiqc_files = ch_reports
    //     .filter { it.exists() }
    // MULTIQC (ch_multiqc_files.collect())
    // multiqc_report = MULTIQC.out.report.toList()
    ch_multiqc_files = Channel
                            .empty()
                            .mix(ch_reports.collect())
    MULTIQC (ch_multiqc_files.collect())
    multiqc_report = MULTIQC.out.report.toList()
}
