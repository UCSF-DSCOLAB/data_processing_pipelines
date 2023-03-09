#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Define the default parameters 
 */
params.input      = "/krummellab/data1/alaa/data/tests/bulk_rnaseq/pipe_dev/test_sample_sheet.csv"   
params.genome     = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
params.genome_idx = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.fai"
params.genome_dict= "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.dict"
params.genome_dir = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/genome_dir"
params.gtf        = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/genomic.gtf" 
params.dbsnp      = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/GCF_000001405.40.gz"
params.dbsnp_tbi  = "/krummellab/data1/alaa/data/human/references/GRCh38.p14/GCF_000001405.40/GCF_000001405.40.gz.tbi"
params.tmp_dir    = "/scratch/alaa/rnax"
params.results    = "results"
// params.gatk       = "/opt/broad/GenomeAnalysisTK.jar"


// Check mandatory parameters (sample sheet)
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' } 
// GATK            =  params.gatk

// Import SUBWORKFLOWS
include { INPUT_CHECK } from './subworkflows/validate_input' 
include { ALIGN_READS } from './subworkflows/align_reads'

// Import MODULES
include { CAT_FASTQ }              from './modules/cat_fastq'
include { GATK4_SPLITNCIGARREADS } from './modules/gatk4_splitncigar'
include { GATK4_BASE_RECALIBRATOR} from './modules/gatk4_recalibrator'
include { GATK4_APPLY_BQSR }       from './modules/gatk4_apply_bqsr'
include { GATK4_HAPLOTYPECALLER }  from './modules/gatk4_haplotype_caller'


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
    // SUBWORKFLOW: Align FastQ reads; sort, and index BAM files
    //
    ALIGN_READS(
        ch_cat_fastq,
        params.gtf,
        params.genome_dir,
        params.tmp_dir
    )
    ch_genome_bam = ALIGN_READS.out.bam
    ch_genome_bai = ALIGN_READS.out.bai
    //
    // MODULE: SplitNCigarReads and reassign mapping qualities
    //
    GATK4_SPLITNCIGARREADS (
        ch_genome_bam,
        ch_genome_bai,
        params.genome,
        params.genome_idx,
        params.genome_dict,
        params.tmp_dir
    )
    ch_split_bam = GATK4_SPLITNCIGARREADS.out.bam
    ch_split_bai = GATK4_SPLITNCIGARREADS.out.bai
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
    ch_bam_bqsr = ch_split_bam.join(ch_recal_table, by: [0])
    ch_bam_bai_bqsr = ch_bam_bqsr.join(ch_split_bai, by: [0])
    GATK4_APPLY_BQSR (
        ch_bam_bai_bqsr,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_bam_variant_calling = GATK4_APPLY_BQSR.out.bam
    ch_bai_variant_calling = GATK4_APPLY_BQSR.out.bai
    ch_bam_bai_variant_calling = ch_bam_variant_calling.join(ch_bai_variant_calling, by: [0])
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_variant_calling,
        params.genome,
        params.genome_idx,
        params.genome_dict
    )
    ch_haplotype_vcf = GATK4_HAPLOTYPECALLER.out.vcf
}