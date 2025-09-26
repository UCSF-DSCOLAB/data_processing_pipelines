# STAR alignment (single rule with static outputs)
rule star_align:
    input:
        reads=lambda wc: reads_for_downstream(wc.sample),
        gtf=GTF,
        genome_dir=GENOME_DIR
    output:
        f"{RESULTS_DIR}/star/{{sample}}.Aligned.sortedByCoord.out.bam",
        f"{RESULTS_DIR}/star/{{sample}}.Aligned.toTranscriptome.out.bam",
        f"{RESULTS_DIR}/star/{{sample}}.ReadsPerGene.out.tab",
        f"{RESULTS_DIR}/star/{{sample}}.Log.final.out"
    threads: 8
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/star
        READS=({input.reads})
        if [ ${{#READS[@]}} -eq 1 ]; then
          READS_CMD="${{READS[0]}}"
        else
          READS_CMD="${{READS[0]}} ${{READS[1]}}"
        fi
        STAR \
            --readFilesIn $READS_CMD \
            --genomeDir {input.genome_dir} \
            --runThreadN {threads} \
            --sjdbGTFfile {input.gtf} \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM GeneCounts \
            --outReadsUnmapped None \
            --outSAMunmapped Within KeepPairs \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:library PL:illumina \
            --outFileNamePrefix {RESULTS_DIR}/star/{wildcards.sample}. \
            --outFilterMismatchNoverLmax {STAR_PARAMS['outfilter_mismatch_n_over_lmax']} \
            --alignSJoverhangMin {STAR_PARAMS['align_sjoverhang_min']} \
            --outFilterMultimapNmax {STAR_PARAMS['outfilter_multimap_nmax']} \
            --seedSearchStartLmax {STAR_PARAMS['seed_search_start_lmax']} \
            {STAR_PARAMS['additional']}
        """