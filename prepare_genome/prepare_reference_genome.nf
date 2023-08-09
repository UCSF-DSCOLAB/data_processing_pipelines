#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Define the default parameters 
*/
params.genome                   = ""
params.gtf                      = ""
params.tmp_dir                  = ""
params.dbsnp                    = ""
params.index_vcf                = false
params.reference_directory      = ""


// Import MODULES
include { SAMTOOLS_FAIDX                     } from './modules/samtools_faidx'
include { STAR_INDEX_GENOME                  } from './modules/star_index_genome'
include { GATK4_GENERATE_SEQUENCE_DICTIONARY } from './modules/gatk4_create_sequence_dictionary'
include { GTF_GENE_FILTER                    } from './modules/gtf_gene_filter'
include { RSEM_PREPAREREFERENCE              } from './modules/rsem_prepare_reference'
include { KALLISTO_INDEX                     } from './modules/kallisto_index'
include { GATK4_INDEX_VCF                    } from './modules/gatk4_index_vcf'


workflow {
    //
    // MODULE: Index the reference genome
    //
    ch_fasta_fai = Channel.empty()
    SAMTOOLS_FAIDX(
        params.genome
    )
    ch_fasta_fai = SAMTOOLS_FAIDX.out.genome_samtools_index
    //
    // MODULE: Create Sequence Dictionary
    //
    ch_gatk_dict = Channel.empty()
    GATK4_GENERATE_SEQUENCE_DICTIONARY (
        params.genome
    )
    ch_gatk_dict = GATK4_GENERATE_SEQUENCE_DICTIONARY.out.genome_dict
    //
    // MODULE: Generate genome indexes using STAR
    //
    ch_star_index = Channel.empty()
    STAR_INDEX_GENOME (
        params.genome,
        params.gtf
    )
    ch_star_index = STAR_INDEX_GENOME.out.index
    //
    // MODULE: Filter GTF file to only include genes
    //
    ch_gtf_genes_only = Channel.empty()
    GTF_GENE_FILTER (
        params.genome,
        params.gtf
    )
    ch_gtf_genes_only = GTF_GENE_FILTER.out.genes_gtf
    //
    // MODULE: Generate genome transcripts using RSEM
    //
    ch_rsem_index = Channel.empty()
    RSEM_PREPAREREFERENCE (
        params.genome,
        ch_gtf_genes_only
    )
    ch_rsem_index = RSEM_PREPAREREFERENCE.out.index
    ch_rsem_transcript = RSEM_PREPAREREFERENCE.out.transcript_fasta
    //
    // MODULE: Generate transcriptome index using Kallisto
    //
    ch_kallisto_index = Channel.empty()
    KALLISTO_INDEX (
        ch_rsem_transcript
    )
    ch_kallisto_index = KALLISTO_INDEX.out.transcript_index
    //
    // MODULE: INDEX VCF IF PROVIDED
    //
    if (params.index_vcf && params.dbsnp) {
        ch_vcf_index = Channel.empty()
        GATK4_INDEX_VCF (
            params.dbsnp
        )
        ch_vcf_index = GATK4_INDEX_VCF.out.tbi
    }
}