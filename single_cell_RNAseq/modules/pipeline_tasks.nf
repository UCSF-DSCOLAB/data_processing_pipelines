

process TEST_GZIP_INTEGRITY {

    input:
    tuple val(library), val(data_type)

    output:
    tuple val(library), val(data_type) 

    """
    gex_library=${library}
    dt=${data_type}
    lib_to_use=\${gex_library/"SCG"/"SC\${dt:0:1}"}

    gzip --test ${params.project_dir}/data/single_cell_${data_type}/raw/\${lib_to_use}/\${lib_to_use}*.fastq.gz
    """
}




/*
 * Step 1a. Run Cellranger
 */
process CELLRANGER {
  
  publishDir "${params.project_dir}/10X_library_csvs/", mode: 'copy', pattern: "${library}_libraries.csv"
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/", pattern: "cellranger/*", mode: 'copy'
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "cellranger.log" }


  container "${params.container.cellranger}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type)
  
  output:
  tuple val(library), path("cellranger/possorted_genome_bam.bam"), path("cellranger/raw_feature_bc_matrix.h5"), emit: bam_h5
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  """

  # create the config
  gex_library=${library}

  echo "fastqs,sample,library_type
    ${params.project_dir}/data/single_cell_GEX/raw/${library},${library},Gene Expression" > ${library}_libraries.csv
  
  # add a line to the config if there is cite-seq data
  if [[ "${data_type}" == "CITE" ]]
  then
      scc_library=\${gex_library/"SCG"/"SCC"}
      echo "${params.project_dir}/data/single_cell_CITE/raw/\${scc_library},\${scc_library},Antibody Capture" >> ${library}_libraries.csv
  fi
  
  cellranger count --id=${library}  \
    --libraries=${library}_libraries.csv \
    --feature-ref=${params.ref.cite_feature_ref} \
    --transcriptome=${params.ref.transcriptome} \
    --localcores=${task.cpus - 1} \
    --localmem=${task.memory.toGiga() - 2}

  mv ${library}/outs cellranger
  """
} 


/*
 * Step 1b. Run Cellranger vdj
 */
process CELLRANGER_VDJ {
  publishDir "${params.project_dir}/data/single_cell_${data_type}/processed/${vdj_library}/", pattern: "cellranger/*", mode: 'copy'
  publishDir "${params.project_dir}/data/single_cell_${data_type}/logs/${vdj_library}/", mode: 'copy', pattern: ".command.log", saveAs: { 
    filename -> "cellranger.log" }

  container "${params.container.cellranger}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type), val(vdj_library) 
  
  output:
  tuple val(library), val(data_type), path("cellranger/clonotypes.csv"), path("cellranger/all_contig_annotations.csv"), emit: vdj_csvs
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  
  """
  vdj_path=${params.project_dir}/data/single_cell_${data_type}/raw/${vdj_library}
  dt=${data_type}
    # TODO: update so that this only occurs on retries if there is a chain error
  if [ \${dt} == "TCR" ]
  then
    chain_type="TR"
  elif [ \${dt} == "BCR" ] 
  then
    chain_type="IG"
  else 
    echo "Warning - chain type should be one of TCR or BCR"
    chain_type="auto"
  fi
  
  cellranger vdj --id=${vdj_library}  \
    --fastqs=\${vdj_path} \
    --reference=${params.ref.vdj_ref} \
    --chain=\${chain_type} \
    --localcores=${task.cpus - 1} \
    --localmem=${task.memory.toGiga() - 2}
  
  mv ${vdj_library}/outs cellranger
  """
}


/*
 * Step 2a. Create barcode list
 */ 
process FILTER_BARCODES{
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "${library}_barcodes_of_interest.filt.list"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
  saveAs: { filename -> "filter_bc.log" }
  
  container "${params.container.rsinglecell}"
  
  input:
  tuple val(library), path(raw_h5) 
  
  output:
  tuple val(library), path("${library}_barcodes_of_interest.filt.list"), emit: bc_list
  path(".command.log"), emit: log
  """
  Rscript ${projectDir}/bin/make_valid_barcodelist.R ${raw_h5} ${params.settings.minfeature} ${params.settings.mincell}
  mv barcodes_of_interest.filt.list ${library}_barcodes_of_interest.filt.list
  
  """
}



/*
 * alternately, filter barcodes based on QC cutoffs
 */
process SEURAT_PRE_FMX_QC {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "pre_fmx_qc.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter", mode: 'copy', pattern: "${library}*"

  container "${params.container.rsinglecell}" 
  containerOptions "-B ${params.settings.default_qc_cuts_dir}"
  
  input:
  tuple val(library), val(data_type), path(raw_h5)

  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library), path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}_diagnostic_plots_pre.pdf"), emit: diagnostic_plots
  path("${library}_quantiles_pre.tsv"), emit: quantiles
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} "null" ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}

  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${data_type} ${library}_initial_raw.rds ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5}
  
  """
}




process SEURAT_PRE_FMX_FILTER {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "pre_fmx_filter.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}_cutoffs.csv"

  container "${params.container.rsinglecell}" 

  input:
  tuple val(library), path(cutoffs), path(raw_sobj)

  output:
  tuple val(library), path("${library}_barcodes_of_interest.filt.list"), emit: bc_list
  tuple val(library), path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}*"), emit: filter_files
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/filter_barcodes_pre_fmx.R ${library} ${projectDir} 
  
  """
}



/*
 * Step 2b. Filter bam file
 */
process FILTER_BAM {
  container "${params.container.popscle}"
  containerOptions "-B ${params.ref.fmx_dir}"

  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "filter_bam.log" }

  input:
  tuple val(library), path(cr_bam), path(barcodes) 

  output:
  tuple val(library), path("${library}_filtered.bam"), emit: bam_file
  path(".command.log"), emit: log
  
  """
  bash filter_bam_for_dsc_pileup.sh \
    ${cr_bam} \
    ${barcodes} \
    ${params.ref.one_k_genome_vcf} \
    ${library}_filtered.bam
  
  """
}

/*
 * Step 2c. Run DSC-pileup
 */
process DSC_PILEUP{
  container "${params.container.popscle}"
  containerOptions "-B ${params.ref.fmx_dir}"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "dsc_pileup.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/${params.settings.demux_method}/", mode: 'copy', pattern: "${library}*.gz"

  input:
  tuple val(library), path(barcodes), path(filtered_bam)

  output:
  tuple val(library), path("${library}*.gz"), emit: plp_files
  path(".command.log"), emit: log

  """
  popscle dsc-pileup --sam ${filtered_bam} \
                   --tag-group CB \
                   --tag-UMI UB \
                   --vcf ${params.ref.one_k_genome_vcf} \
                   --group-list ${barcodes} \
                   --out ${library}
  """
}


/*
 * Step 2d. Merge DSC pileups
 */
process MERGE_DSC { 
  publishDir "${params.settings.merge_demux_dir}/${pool}/", mode: 'copy', pattern: "${pool}.tsv"
  publishDir "${params.settings.merge_demux_dir}/${pool}", mode: 'copy', pattern: "${pool}*.gz"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${pool}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "merge_fmx.log" }

  container "${params.container.python}"

  input:
  // TODO: expand out what these path files are...
  tuple val(pool), path(pool_files)
  
  output:
  tuple val(pool), path("${pool}.tsv"), path("merged.plp.gz"), path("merged.var.gz"), path("merged.cel.gz"), path("${pool}.barcodes.gz"), emit: merged_files
  path(".command.log"), emit: log

  """
  source utils.sh

  # write out the config file
  printf "sample\n" > ${pool}.tsv
  echo `ls *.plp.gz`
  for my_dsc in `ls *.plp.gz`
  do
    library=\$(basename \${my_dsc%%.*})
    printf "\${library}\n" >> ${pool}.tsv
  done
  
  # merge
  python ${projectDir}/bin/merge_freemuxlet_dsc_pileups.py --freemuxlet_dir_tsv ${pool}.tsv --ignore_diff_lengths_error
  gzipStripDttm merged.plp.gz
  gzipStripDttm merged.cel.gz
  gzipStripDttm merged.var.gz
  gzipStripDttm merged.barcodes.gz && mv merged.barcodes.gz ${pool}.barcodes.gz
  """
}

/*
 * Step 2e. Run Freemuxlet for a pool
 */ 
process FREEMUXLET_POOL {
  publishDir "${params.settings.merge_demux_dir}/${pool}/", 
    mode: 'copy', 
    pattern: "merged*", 
    saveAs : { filename -> "${pool}_${filename}" }

  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${pool}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "freemuxlet.log" }

  container "${params.container.popscle}"

  input:
  tuple val(pool), val(nsamples), path(merged_plp), path(merged_var), path(merged_cel), path(merged_barcodes)
  
  output:
  tuple val(pool), path("merged.clust1.samples.gz"), path("merged.lmix"), path('merged.clust1.vcf.gz'), emit: merged_files
  tuple val(pool), path("merged.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """
  source utils.sh

  # run freemuxlet
  popscle freemuxlet --plp merged \
                     --out merged \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed} \
                     --group-list ${merged_barcodes}
  
  gzipStripDttm merged.clust1.samples.gz
  vcfStripDttm merged.clust1.vcf.gz
  
  """
} 

/*
 * Step 2e. Run Freemuxlet for a library
 */ 
process FREEMUXLET_LIBRARY {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "freemuxlet.log" }

  container "${params.container.popscle}"

  input:
  tuple val(library), val(nsamples), path(plp_files)
  
  output:
  tuple val(library), path("${library}*"), emit: fmx_files
  tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
  tuple val(library), path("${library}.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """
  source utils.sh

  popscle freemuxlet --plp ${library} \
                     --out ${library} \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed} 
  
  # then unzip and do next steps   
  vcfStripDttm ${library}.clust1.vcf.gz
  gzipStripDttm ${library}.plp.gz
  gzipStripDttm ${library}.cel.gz
  gzipStripDttm ${library}.var.gz
  gzipStripDttm ${library}.umi.gz

  demuxTsvFromFmx ${library} 
  """
} 


/*
 * Run assign to gt for a pool or library 
 */
process FMX_ASSIGN_TO_GT {
  publishDir "${params.project_dir}/fmx_assign_to_gt/${pool}/", mode: 'copy', pattern: "${pool}_gtcheck*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${pool}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "fmx_assign_to_gt.log" }

  container "${params.container.rplus_bcftools}"
  containerOptions "-B ${params.project_dir} -B ${params.settings.ref_vcf_dir}"

  input:
  tuple val(pool), val(ref_vcf), path(fmx_vcf) 

  output:
  tuple val(pool), path("${pool}*"), emit: outfiles
  path(".command.log"), emit: log


  """
  bash ${projectDir}/bin/run_gtcheck.sh ${pool} ${params.settings.ref_vcf_dir}/${ref_vcf} ${fmx_vcf} ${params.settings.ref_vcf_type}
  Rscript ${projectDir}/bin/examine_gtcheck.R ${pool} ${pool}_gtcheck.out
  """

}


/*
 * Run demuxlet for a library
 */
process DEMUXLET_LIBRARY {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/demuxlet", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "demuxlet.log" }
  container "${params.container.popscle}"
  
  input:
  tuple val(library), path(vcf), path(plp_files)

  output:
  tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
  path("${library}.clust1.samples.gz"), emit: dmx_files
  path(".command.log"), emit: log
  """
  source utils.sh

  popscle demuxlet --plp ${library} \
                     --out ${library} \
                     --vcf ${vcf} \
                     --field GT 
  
  # select desired files              
  extractDemuxTsv ${library}.best ${library}
  mv ${library}.best ${library}.clust1.samples && gzip -f -n ${library}.clust1.samples
  """
}

process DEMUXLET_POOL {
 publishDir "${params.settings.merge_demux_dir}/${pool}", 
    mode: 'copy', 
    pattern: "merged*", 
    saveAs : { filename -> "${pool}_${filename}" }

  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${pool}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "demuxlet.log" }
  container "${params.container.popscle}"
  
  input:
  tuple val(pool), path(vcf), path(merged_plp), path(merged_var), path(merged_cel), path(merged_barcodes)

  output:
  tuple val(pool), path("merged.best"), emit: merged_best
  path(".command.log"), emit: log

  """
  popscle demuxlet --plp merged \
                     --out merged \
                     --vcf ${vcf} \
                     --field GT 
  
   """
}

process SEPARATE_DMX {
   publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/demuxlet", mode: 'copy', pattern: "${library}*"
   publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
   saveAs: { filename -> "separate_dmx.log" }


  input: 
   tuple val(library), path(merged_best)
  output:
   tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path("${library}.clust1.samples.gz"), emit: full_sample_file
   path(".command.log"), emit: log

  """
  source utils.sh

  head -1 ${merged_best} > ${library}.clust1.samples0
  grep "${library}\\s" ${merged_best} >> ${library}.clust1.samples0
  sed \"s/--${library}//g\" ${library}.clust1.samples0 > ${library}.clust1.samples
  extractDemuxTsv ${library}.clust1.samples ${library}
  gzip -f -n ${library}.clust1.samples
  """
}





process UNMERGE_FMX {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${pool}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "unmerge_fmx.log" }

  container "${params.container.python}"

  input:
  tuple val(pool), path(fmx_tsv), path(merged_samples), path(merged_lmix), path(merged_vcf)
  
  output:
  tuple path("*vcf.gz"), path("*samples.gz"), path("*.lmix"), emit: samples_file
  path(".command.log"), emit: log

  """
  # unmerge
  python ${projectDir}/bin/unmerge_freemuxlet_dsc_pileups.py --freemuxlet_dir_tsv ${fmx_tsv} 
  rm merged*
  """
}

process SEPARATE_FMX {
   publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "${library}*"
   publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "separate_fmx.log" }
  input:
   tuple val(library), path(vcf_file), path(sample_file), path(lmix_file)

  output:
   tuple val(library), path("${library}.clust1.samples.gz"), path("${library}.clust1.vcf.gz"), path("${library}.lmix"), emit: fmx_files
   tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path(".command.log"), emit: log

  """
  source utils.sh

  demuxTsvFromFmx ${library}
  vcfStripDttm ${library}.clust1.vcf.gz

  """
}


/*
 * Step 3. Run DoubletFinder
 */
process FIND_DOUBLETS {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/finding_doublets", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "run_df.log" }


  container "${params.container.rsinglecell}" 

  input:
  tuple val(library), val(ncells_loaded), path(raw_h5), path(fmx_clusters)

  output:
  tuple val(library), path("${library}_seurat_object_findingDoublets.rds"), emit: sobj
  path("${library}*"), emit: df_stats
  path(".command.log"), emit: log

  """
  
  Rscript ${projectDir}/bin/find_doublets.R ${raw_h5} ${fmx_clusters} ${library} ${ncells_loaded} ${params.settings.minfeature} ${params.settings.mincell} ${params.settings.randomseed} ${projectDir}
  
  """
}

process LOAD_SOBJ {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "load_sobj.log" }

  container "${params.container.rsinglecell}" // todo update to one with DF

  input:
  tuple val(library), path(raw_h5), path(fmx_clusters)

  output:
  tuple val(library), path("${library}*.rds"), emit: sobj
  path(".command.log"), emit: log

  """
  
  Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} ${fmx_clusters} ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}
  
  """
}








/* 
 * Step 5a. Seurat add vdj
 */
process SEURAT_ADD_BCR {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "add_bcr.log" }

  container "${params.container.rsinglecell}"

  
  input:
  tuple val(library), path(sobj), val(data_type), path(clonotypes_csv), path(contig_csv)
  
  output:
  tuple val(library), path("${library}_w_BCR.RDS"), emit: sobj
  path(".command.log"), emit: log

  """
  if [[ "${data_type}" == "no BCR" ]]
  then
    cp ${sobj} "${library}_w_BCR.RDS" # todo - switch to soft link
  else
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${clonotypes_csv} ${contig_csv} ${projectDir}
  fi
  
  """

}

/* 
 * Step 5b. Seurat add vdj
 */
process SEURAT_ADD_TCR {
    
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "add_tcr.log" }

  container "${params.container.rsinglecell}"

  
  input:
  tuple val(library), path(sobj), val(data_type), path(clonotypes_csv), path(contig_csv)
  
  output:
  tuple val(library), path("${library}_w_TCR.RDS"), emit: sobj
  path(".command.log"), emit: log

  """
  if [[ "${data_type}" == "no TCR" ]]
  then
    cp ${sobj} "${library}_w_TCR.RDS" # todo - switch to soft link
  else
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${clonotypes_csv} ${contig_csv} ${projectDir}
  fi
  
  """

}

/* 
 * Step 6. Run Seurat
 */
process SEURAT_QC {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "seurat_qc.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"

  container "${params.container.rsinglecell}"
  containerOptions "-B ${params.settings.default_qc_cuts_dir}"

  input:
  tuple val(library), val(main_dt), path(doublet_finder_sobj), path(raw_h5)
  
  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}_diagnostic_plots_pre.pdf"), emit: diagnostic_plots
  path("${library}_quantiles_pre.tsv"), emit: quantiles
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5}

  """
}

/* 
 * Step 6b. Run Seurat
 */
process SEURAT_LOAD_POST_QC {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "seurat_qc.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"
  // For testing
  publishDir "${workDir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"

  container "${params.container.rsinglecell}"

  input:
  tuple val(library), val(main_dt), path(doublet_finder_sobj), path(raw_h5), path(cutoffs)
  
  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${cutoffs} ${raw_h5}

  """
}


/* 
 * Step 5b. Run Post filter
 */
process SEURAT_POST_FILTER {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "post_filter.log" }
    
  container "${params.container.rsinglecell}"

  input:
  tuple val(library), path(cutoffs), path(raw_sobj), path(raw_h5)
  
  output:
  tuple val(library), path("${library}_filtered.rds"), emit: post_process
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}*"), emit: outfiles
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/process_with_seurat_post_filter.R ${library} ${projectDir} ${raw_h5} ${params.settings.remove_demux_DBL} ${params.settings.remove_all_DBL}
  """
}

