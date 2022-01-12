#!/usr/bin/env nextflow


// This Pipeline expects the following to be provided in yml form on the command line
/* params.samples is a map of
 *       sample_name:
 *           - XXX_R1_YYY.fastq.gz
 *           - XXX_R2_YYY.fastq.gz
 */
// params.outdir is a string value
// params.cohort_name is a string value(e.g. plate31)

log.info """\
         DSCo/Genomics RNASeq Pipeline
         =============================
         genome: ${params.samples}
         annot : ${params.outdir}
         reads : ${params.cohort_name}
         """
         .stripIndent()


// https://github.com/nf-core/vipr/blob/master/main.nf
def GetReadPair = { sk ->
    tuple(file(params.samples[sk]['R1']),
          file(params.samples[sk]['R2']))
}

/* FIXME allow other means of defining input, e.g. CSV.
 * input for fastq files. channel has sample name as key and all read pairs following
 * see https://groups.google.com/forum/#!topic/nextflow/CF7Joh5xrkU
 */
sample_keys = params.samples.keySet()
println "List of samples: " +  sample_keys.join(", ")
Channel
    .from(sample_keys)
    .map { sk -> tuple(sk, GetReadPair(sk)) }
    .set { readMapChannel }

/*
 * Step 1. Test Gzip Integrity
 */
process testGzipIntegrity {
    tag { "${sample}--${params.cohort_name}" }

    input:
    tuple sample, path(reads) from readMapChannel

    output:
    tuple sample, path(reads) into gzipOKReads

    """
    gzip --test ${reads[0]} ${reads[1]}
    """
}

/*
 * Step 2. Calculate md5sum
 */
process md5SumCalculation {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/qc/", mode: 'copy', pattern: "${sample}_fastq_md5sums.txt"

    input:
    tuple sample, path(reads) from gzipOKReads

    output:
    tuple sample, path(reads) into md5dReads
    path "${sample}_fastq_md5sums.txt"

    """
    md5sum ${reads[0]} ${reads[1]} > ${sample}_fastq_md5sums.txt
    """
}

/*
 * Step 3. Run fastp to trim adapters
 */
process runFastp {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/qc/", mode: 'copy', pattern: "${sample}.fastp.*"

    container "${params.container.fastp}"
    cpus 12
    memory '10G'

    input:
    tuple val(sample), path(reads) from md5dReads

    output:
    tuple val(sample), path("${sample}_trimmed_R{1,2}.fq.gz") into trimmedReads
    path "${sample}.fastp.json"
    path "${sample}.fastp.html"

    """
    fastp -i ${reads[0]} \
          -o ${sample}_trimmed_R1.fq.gz \
          -I ${reads[1]} \
          -O ${sample}_trimmed_R2.fq.gz \
          --length_required 20 \
          --adapter_sequence ${params.ref.adapter_seq1} \
          --adapter_sequence_r2 ${params.ref.adapter_seq2} \
          --correction  \
          --trim_poly_g  \
          --thread ${task.cpus} \
          -j ${sample}.fastp.json \
          -h ${sample}.fastp.html
    """
}



/*
 * Step 4. Align to the human genome/transcriptome
 */
process alignToGenomeTranscriptome {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Chimeric.out.junction"
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.ReadsPerGene.out.tab"
    
    container "${params.container.star}"
    containerOptions "-B ${params.ref.star_genome_dir}"
    cpus 32
    memory '100G'

    input:
    tuple val(sample), path(reads) from trimmedReads

    output:
    tuple val(sample), "${sample}.trimmed.star.Unmapped.out.mate?" into starUnamppedReads
    tuple val(sample), "${sample}.trimmed.star.Aligned.toTranscriptome.out.bam" into transcriptomeBAM
    tuple val(sample), "${sample}.trimmed.star.Aligned.sortedByCoord.out.bam" into genomeBAM
    path "${sample}.trimmed.star.Chimeric.out.junction"
    path "${sample}.trimmed.star.ReadsPerGene.out.tab"

    """
    STAR --readFilesIn ${reads[0]} ${reads[1]} \
         --genomeDir ${params.ref.star_genome_dir} \
         --outSAMunmapped None \
         --outFilterType BySJout \
         --outSAMattributes NH HI AS NM MD \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMax 100000 \
         --alignMatesGapMax 100000 \
         --alignSJDBoverhangMin 10 \
         --sjdbScore 1 \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMheaderHD @HD VN:1.4 SO:coordinate \
         --readFilesCommand zcat \
         --alignSJstitchMismatchNmax 5 -1 5 5 \
         --chimSegmentMin 12 \
         --chimJunctionOverhangMin 12 \
         --chimSegmentReadGapMax 3 \
         --chimMultimapScoreRange 10 \
         --chimMultimapNmax 10 \
         --chimNonchimScoreDropMin 10 \
         --outReadsUnmapped Fastx \
         --outSAMattrRGline ID:GRPundef \
         --outSAMstrandField intronMotif \
         --peOverlapNbasesMin 12 \
         --peOverlapMMp 0.1 \
         --runThreadN ${task.cpus} \
         --twopassMode Basic \
         --outFileNamePrefix ${sample}.trimmed.star.
    
    """
}



/*
 * Step 5. gzip the STAR unmapped reads for storage
 */
process gzipSTARUnmappedReads {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Unmapped.out.mate?.fq.gz"

    container "${params.container.utils}"
    cpus 16
    memory '50G'

    input:
    tuple val(sample), path(reads) from starUnamppedReads

    output:
    path "${sample}.trimmed.star.Unmapped.out.mate1.fq.gz"
    path "${sample}.trimmed.star.Unmapped.out.mate2.fq.gz"

    """
    pigz -cf -p ${task.cpus} ${reads[0]} > ${sample}.trimmed.star.Unmapped.out.mate1.fq.gz
    pigz -cf -p ${task.cpus} ${reads[1]} > ${sample}.trimmed.star.Unmapped.out.mate2.fq.gz
    """
}

/*
 * Step 6. Convert STAR transcriptome bam for storage
 */
process transcriptomeBAMToCRAM {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.cram"

    container "${params.container.samtools}"
    containerOptions "-B ${params.ref.rsem_star_transcript_ref_dir}"
    cpus 16
    memory '50G'

    input:
    tuple val(sample), path(t_bam) from transcriptomeBAM

    output:
    tuple val(sample), path(t_bam) into transcriptomeBAM2
    path "${sample}.trimmed.star.cram"

    """
    samtools view  -@ ${task.cpus} \
                   -C \
                   --no-PG \
                   -T ${params.ref.rsem_star_transcript_ref_dir}/${params.ref.rsem_star_transcript_ref_prefix}.transcripts.fa \
                   -o ${sample}.trimmed.star.cram \
                   ${t_bam}
    """
}

/*
 * Step 7. MarkDuplicates on the STAR genome Bam
 */
process markDuplicatesGenomicBAM {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Aligned.sortedByCoord.out.duplication_metrics"

    container "${params.container.picard}"
    cpus 1
    memory '50G'

    input:
    tuple val(sample), path(g_bam) from genomeBAM

    output:
    tuple val(sample), "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.bam" into dedupGenomeBAM
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.duplication_metrics"

    """
    java -Xmx20g \
         -jar /opt/picard/picard.jar \
                MarkDuplicates \
                -VALIDATION_STRINGENCY SILENT \
                -REMOVE_DUPLICATES false \
                -ASSUME_SORTED true \
                -INPUT ${g_bam} \
                -OUTPUT ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.bam \
                -METRICS_FILE ${sample}.trimmed.star.Aligned.sortedByCoord.out.duplication_metrics
    """
}

/*
 * Step 8. Convert Deduplicated STAR genome bam for storage
 */
process dedupGenomeBAMToCRAM {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram*"

    container "${params.container.samtools}"
    containerOptions "-B ${params.ref.sequence_ref_dir}"
    cpus 16
    memory '50G'

    input:
    tuple val(sample), path(dg_bam) from dedupGenomeBAM

    output:
    tuple val(sample), path(dg_bam) into dedupGenomeBAM2
    tuple val(sample), path(dg_bam) into dedupGenomeBAM3
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram"
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram.crai"
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram.flagstat"

    """
    samtools view  -@ ${task.cpus} \
                   --no-PG \
                   --write-index \
                   -C \
                   -T ${params.ref.sequence_ref_dir}/${params.ref.fasta_file} \
                   -o ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram \
                   ${dg_bam}

    samtools flagstat  ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram \
        > ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.cram.flagstat
    """
}

/*
 * Step 9. Quantify Transcripts using RSEM
 */
process runRSEM {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.rsem.*.results"

    container "${params.container.rsem}"
    containerOptions "-B ${params.ref.rsem_star_transcript_ref_dir}"
    cpus 12
    memory '50G'

    input:
    tuple val(sample), path(t_bam) from transcriptomeBAM2

    output:
    path "${sample}.rsem.genes.results"
    path "${sample}.rsem.isoforms.results"
    
    """
    rsem-calculate-expression --paired-end  \
                              --no-bam-output  \
                              --alignments  \
                              --num-threads ${task.cpus} \
                              ${t_bam} \
                              ${params.ref.rsem_star_transcript_ref_dir}/${params.ref.rsem_star_transcript_ref_prefix} \
                              ${sample}.rsem   
    """
}

/*
 * Step 10. Get CollectRnaSeqMetrics on the dedup genome BAM
 */
process dedupGenomeBAMRSQMetrics {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.rnaseq_metrics"

    container "${params.container.picard}"
    containerOptions "-B ${params.ref.genome_flat_ref_dir} -B ${params.ref.ribosomal_intervals_dir}"
    cpus 12
    memory '50G'

    input:
    tuple val(sample), path(dg_bam) from dedupGenomeBAM2

    output:
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.rnaseq_metrics"
    
    """
    java -Xmx${task.memory.toGiga()-5}g \
         -jar /opt/picard/picard.jar \
         CollectRnaSeqMetrics \
            -VALIDATION_STRINGENCY SILENT \
            -ASSUME_SORTED true \
            -STRAND_SPECIFICITY NONE \
            -REF_FLAT ${params.ref.genome_flat_ref_dir}/${params.ref.genome_flat_ref_filename}  \
            -RIBOSOMAL_INTERVALS ${params.ref.ribosomal_intervals_dir}/${params.ref.ribosomal_intervals_filename}  \
            -I ${dg_bam} \
            -OUTPUT ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.rnaseq_metrics
    """
}

/*
 * Step 11. Get CollectAlignmentSummaryMetrics on the dedup genome BAM
 */
process dedupGenomeBAMAlignmentMetrics {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.alignment_metrics"

    container "${params.container.picard}"
    cpus 12
    memory '50G'

    input:
    tuple val(sample), path(dg_bam) from dedupGenomeBAM3

    output:
    path "${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.alignment_metrics"

    """
    java -Xmx${task.memory.toGiga()-5}g \
         -jar /opt/picard/picard.jar \
         CollectAlignmentSummaryMetrics \
            -VALIDATION_STRINGENCY SILENT \
            -INPUT ${dg_bam} \
            -OUTPUT ${sample}.trimmed.star.Aligned.sortedByCoord.out.deduplicated.alignment_metrics
    """
}


workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}