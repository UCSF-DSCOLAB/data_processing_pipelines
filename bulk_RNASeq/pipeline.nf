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
    tuple val(sample), path("${sample}_trimmed_R{1,2}.fq.gz") into trimmedReads2
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
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.star.Aligned.toTranscriptome.out.bam"    

    container "${params.container.star}"
    containerOptions "-B ${params.ref.star_genome_dir}"
    cpus 32
    memory '100G'

    input:
    tuple val(sample), path(reads) from trimmedReads

    output:
    tuple val(sample), "${sample}.trimmed.star.Aligned.toTranscriptome.out.bam" into transcriptomeBAM
    tuple val(sample), "${sample}.trimmed.star.Aligned.sortedByCoord.out.bam" into genomeBAM
    tuple val(sample), "${sample}.trimmed.star.Aligned.sortedByCoord.out.bam" into genomeBAM2
    path "${sample}.trimmed.star.Chimeric.out.junction"
    path "${sample}.trimmed.star.ReadsPerGene.out.tab"

    """
    STAR --readFilesIn ${reads[0]} ${reads[1]} \
         --genomeDir ${params.ref.star_genome_dir} \
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
         --outReadsUnmapped None \
	 --outSAMunmapped Within KeepPairs \
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
 * Step 5. Extract unmapped reads and gzip them
 */
process extractUnmappedReads {

    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.trimmed.unmapped.{1,2}.fq.gz" 

    container "${params.container.picard}"
    cpus 16
    memory '50G'

    input: 
     tuple val(sample), path(bamfile) from genomeBAM
    
    output:
     path "${sample}.trimmed.unmapped.{1,2}.fq.gz" 

    """
    java -Xmx${task.memory.toGiga()-5}g \
         -jar /opt/picard/picard.jar \
            ViewSam \
                VALIDATION_STRINGENCY=SILENT \
                ALIGNMENT_STATUS=Unaligned \
                PF_STATUS=All \
                I=${bamfile} \
                > ${sample}.trimmed.unmapped.bam

    java -Xmx${task.memory.toGiga()-5}g \
         -jar /opt/picard/picard.jar \
            SamToFastq \
                VALIDATION_STRINGENCY=SILENT \
                I=${sample}.trimmed.unmapped.bam \
                FASTQ=${sample}.trimmed.unmapped.1.fq.gz \
                SECOND_END_FASTQ=${sample}.trimmed.unmapped.2.fq.gz

    """
}


/*
 * Step 6. Extract and convert mapped transcriptome
 */
process transcriptomeBAMToCRAM {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.transcriptome.trimmed.star.cram"

    container "${params.container.samtools}"
    containerOptions "-B ${params.ref.rsem_star_transcript_ref_dir}"
    cpus 16
    memory '50G'

    input:
    tuple val(sample), path(t_bam) from transcriptomeBAM

    output:
    tuple val(sample), "${sample}.trimmed.star.Transcriptome.mapped.bam" into transcriptomeBAM2
    path "${sample}.trimmed.star.Transcriptome.mapped.cram"


    """
    samtools view  -b \
		   -F 0x4 \
    	     	   -@ ${task.cpus} \
                   --no-PG \
		   -o ${sample}.trimmed.star.Transcriptome.mapped.bam \
		   ${t_bam}

    samtools view -@ ${task.cpus} \
		   -C \
		   --no-PG \
                   -T ${params.ref.rsem_star_transcript_ref_dir}/${params.ref.rsem_star_transcript_ref_prefix}.transcripts.fa \
                   -o ${sample}.trimmed.star.Transcriptome.mapped.cram \
                   ${sample}.trimmed.star.Transcriptome.mapped.bam
    """
}


/*
 * Step 7b. Extract and convert mapped genome
 */
process extractMappedGenome {
    tag { "${sample}--${params.cohort_name}" }


    container "${params.container.samtools}"
    containerOptions "-B ${params.ref.rsem_star_transcript_ref_dir}"
    cpus 16
    memory '50G'

    input:
    tuple val(sample), path(g_bam) from genomeBAM2

    output:
    tuple val(sample), "${sample}.trimmed.star.Aligned.mapped.bam" into genomeBAM3


    """
    samtools view  -b \
		   -F 0x4 \
    	     	   -@ ${task.cpus} \
                   --no-PG \
		   -o ${sample}.trimmed.star.Aligned.mapped.bam \
		   ${g_bam}

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
    tuple val(sample), path(g_bam) from genomeBAM3

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
    tuple val(sample), path(dg_bam) into dedupGenomeBAM4
    tuple val(sample), path(dg_bam) into dedupGenomeBAM5
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
 * Step 10. Use kallisto to quantify reads
 */
process runKallisto {
    tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.kallisto.abundance.*"
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.kallisto.run_info.json"

    container "${params.container.kallisto}"
    containerOptions "-B ${params.ref.kallisto_dir}"
    cpus 12
    memory '50G'

    input:
    tuple val(sample), path(reads) from trimmedReads2

    output:
    path "${sample}.kallisto.abundance.h5"
    path "${sample}.kallisto.abundance.tsv"
    path "${sample}.kallisto.run_info.json"
    
    """
    kallisto quant -i ${params.ref.kallisto_dir}/${params.ref.kallisto_index} -o ${sample} -t ${task.cpus} ${reads[0]} ${reads[1]}
    mv ${sample}/abundance.h5 ${sample}.kallisto.abundance.h5
    mv ${sample}/abundance.tsv ${sample}.kallisto.abundance.tsv
    mv ${sample}/run_info.json ${sample}.kallisto.run_info.json

    """
}


/*
 * Step 11. Get CollectRnaSeqMetrics on the dedup genome BAM
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
 * Step 12. Get CollectAlignmentSummaryMetrics on the dedup genome BAM
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


/*
 * Step 13a: add read groups
 */
process addRGroups {
  tag { "${sample}--${params.cohort_name}" }
 
  container "${params.container.picard}"
  containerOptions "-B ${params.dirs.genome}"

  cpus 16
  memory '50G'

  input: tuple val(sample), path(bamfile) from dedupGenomeBAM4

  output: tuple val(sample), "${sample}.withRG.bam" into rgBAM

  """
  java -Xmx${task.memory.toGiga()-5}g \
              -jar /opt/picard/picard.jar \
	      AddOrReplaceReadGroups \
                                        VALIDATION_STRINGENCY=SILENT \
                                        CREATE_INDEX=TRUE \
                                        INPUT=${bamfile} \
                                        OUTPUT=${sample}.withRG.bam \
                                        SORT_ORDER=coordinate \
                                        RGID=grp1 \
                                        RGPL="ILLUMINA" \
                                        RGSM=grp \
                                        RGLB=grpA \
                                        RGPU=grpB
 					
  """
}

/*
 * Step 13c: filter/split BAM pre-GATK
 */
process splitBAM {
   tag { "${sample}--${params.cohort_name}" }

   container "${params.container.gatk}"
   containerOptions "-B ${params.dirs.genome}"
   
   cpus 16
   memory '50G'

   input: 
     tuple val(sample), path(rg_bam) from rgBAM
   
   output: 
     tuple val(sample), "${sample}.withRG.nSplit.bam" into splitBAM


   """
   gatk SplitNCigarReads \
     -R ${params.ref.sequence_ref_dir}/${params.ref.fasta_file} \
     -I ${rg_bam} \
     -O ${sample}.withRG.nSplit.bam
	
   """
}

/*
 * Step 14: run GATK
 */
process runGATK {
    tag { "${sample}--${params.cohort_name}" }

   container "${params.container.gatk}"
   containerOptions "-B ${params.dirs.genome}"

   cpus 16
   memory '50G'

   input: 
     tuple val(sample), path(split_bam) from splitBAM
   output: 
     tuple val(sample), "${sample}.raw.vcf" into gatkOut
     tuple val(sample), "${sample}.raw.vcf.idx" into gatkOutIdx
  """
  gatk HaplotypeCaller \
    -L ${params.ref.exon_bed} \
    -R ${params.ref.sequence_ref_dir}/${params.ref.fasta_file} \
    -D ${params.ref.dbsnp} \
    -I ${split_bam} \
    --dont-use-soft-clipped-bases true -stand-call-conf 20.0 \
    -O ${sample}.raw.vcf
  """

}

/*
 * Step 15: filter based on GATK
 */
process filterGATK {
   tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.filtered.vcf"    
  container "${params.container.gatk}"
  containerOptions "-B ${params.dirs.genome}"

  cpus 16
  memory '50G'

  input: 
     tuple val(sample), path(raw_gatk) from gatkOut
     tuple val(sample), path(raw_gatk_idx) from gatkOutIdx
  output: 
    path "${sample}.filtered.vcf"

 """
 gatk VariantFiltration \
    -L ${params.ref.exon_bed} \
    -R ${params.ref.sequence_ref_dir}/${params.ref.fasta_file} \
    -V ${raw_gatk} \
    -window 35 -cluster 3 \
    --filter-name FS -filter "FS > 30.0" \
    --filter-name QD -filter "QD < 2.0" \
    -O ${sample}.filtered.vcf
  """
}

/*
 * Step 16: Run arcasHLA
 */
process arcasHLA {
   tag { "${sample}--${params.cohort_name}" }
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.genes.json"
    publishDir "${params.outdir}/${sample}/alignments/", mode: 'copy', pattern: "${sample}.genotype.json"

 cpus 8
 memory '30G'

 container "${params.container.arcashla}"

 input: tuple val(sample), path(bamfile) from dedupGenomeBAM5

 output:
 path "${sample}.genes.json"
 path "${sample}.genotype.json"

 """
 cp "${bamfile}" "${sample}.bam"
 
 arcasHLA extract -t 8 -v "${sample}.bam"
    
 arcasHLA genotype --genes A,B,C,DPB1,DQB1,DQA1,DRB1 \
                        --threads 8 \
                        --verbose \
			"${sample}.extracted.1.fq.gz" \
			"${sample}.extracted.2.fq.gz"
 """

}


workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}