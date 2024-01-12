

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
    
    if [[ "${data_type}" == "CITE" ]]
    then
      gzip --test ${params.project_dir}/data/single_cell_GEX/raw/${library}/${library}*.fastq.gz
    fi

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
  my_dir=\${PWD}
  cd \${TMPDIR}
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
    --transcriptome=${params.ref.transcriptome} 

  mv ${library}/outs \${my_dir}/cellranger
  """
} 


/*
 * Step 1b. Run Cellranger vdj
 */
process CELLRANGER_VDJ {
  publishDir "${params.project_dir}/data/single_cell_${data_type}/processed/${library}/cellranger", mode: 'copy'
  publishDir "${params.project_dir}/data/single_cell_${data_type}/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "cellranger.log" }

  container "${params.container.cellranger}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type) 
  
  output:
  tuple val(library), val(data_type), path("clonotypes.csv"), path("filtered_contig_annotations.csv"), emit: bam_h5
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  
  """
  my_dir=\${PWD}
  cd \${TMPDIR}
  echo \${TMPDIR}

  gex_library=${library}
  data_type_name = "${data_type}"
  vdj_library=\${gex_library/"SCG"/"SC\${data_type_name:0:1}"}

  vdj_path=\${params.project_dir}/data/single_cell_${data_type}/raw/\${vdj_library}
  
  cellranger vdj --id="\${vdj_library}"  \
    --fastqs=\${vdj_path} \
    --reference="\${vdj_library}" \
    --reference=${params.ref.vdj_ref} 
  

  mv ${vdj_library}/outs \${my_dir}/cellranger

  """
}


/*
 * Step 2a. Create barcode list
 */ 
process FILTER_BARCODES{
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "barcodes_of_interest.filt.list"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "filter_bc.log" }
  
  container "${params.container.rsinglecell}"
  
  input:
  tuple val(library), path(raw_h5) 
  
  output:
  tuple val(library), path("barcodes_of_interest.filt.list"), emit: bc_list
  path(".command.log"), emit: log
  """
  Rscript ${projectDir}/bin/make_valid_barcodelist.R ${raw_h5} ${params.settings.minfeature} ${params.settings.mincell}
  
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
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter", mode: 'copy', pattern: "barcodes_of_interest.filt.list"

  container "${params.container.rsinglecell}" 

  input:
  tuple val(library), path(cutoffs), path(raw_sobj)

  output:
  tuple val(library), path("barcodes_of_interest.filt.list"), emit: bc_list
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
  mv merged.barcodes.gz ${pool}.barcodes.gz

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
  containerOptions "-B ${params.ref.fmx_dir}"
  
  input:
  tuple val(pool), val(nsamples), path(merged_plp), path(merged_var), path(merged_cel), path(merged_barcodes)
  
  output:
  tuple val(pool), path("merged.clust1.samples.gz"), path("merged.lmix"), path('merged.clust1.vcf.gz'), emit: merged_files
  tuple val(pool), path("merged.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """
  # run freemuxlet
  popscle freemuxlet --plp merged \
                     --out merged \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed} \
                     --group-list ${merged_barcodes}
  
  
  """
} 

/*
 * Step 2e. Run Freemuxlet for a library
 */ 
process FREEMUXLET_LIBRARY {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "freemuxlet.log" }

  container "${params.container.popscle}"
  containerOptions "-B ${params.ref.fmx_dir}"
  
  input:
  tuple val(library), val(nsamples), path(plp_files)
  
  output:
  tuple val(library), path("${library}*"), emit: fmx_files
  tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
  tuple val(library), path("${library}.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """

  popscle freemuxlet --plp ${library} \
                     --out ${library} \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed} 
  
  # then unzip and do next steps                   
  gunzip -f ${library}.clust1.samples.gz
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} ${library}.clust1.samples > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples

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
  path("${library}*"), emit: dmx_files
  path(".command.log"), emit: log
  """
  popscle demuxlet --plp ${library} \
                     --out ${library} \
                     --vcf ${vcf} \
                     --field GT 
  
  # select desired files              
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} ${library}.best > ${library}.clust1.samples.reduced.tsv
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
  head -1 ${merged_best} > ${library}.clust1.samples0
  grep "${library}\\s" ${merged_best} >> ${library}.clust1.samples0
  sed \"s/--${library}//g\" ${library}.clust1.samples0 > ${library}.clust1.samples
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} "${library}.clust1.samples" > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples
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
   tuple val(library), path(library_files)

  output:
   tuple path("${library}.clust1.samples.gz"), path("${library}.clust1.vcf.gz"), path("${library}.lmix"), emit: fmx_files
   tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path(".command.log"), emit: log

  """
  gunzip -f ${library}.clust1.samples.gz
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} ${library}.clust1.samples > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples
  """
}

// TODO: unify the two processes
process SEPARATE_FMX_PRE {
   publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "${library}*"
   publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "separate_fmx.log" }

  input:
   tuple val(library), path(vcf_file), path(sample_file), path(lmix_file)

  output:
   tuple path("${library}.clust1.samples.gz"), path("${library}.clust1.vcf.gz"), path("${library}.lmix"), emit: fmx_files
   tuple val(library), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path(".command.log"), emit: log

  """
  gunzip -f ${library}.clust1.samples.gz
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} ${library}.clust1.samples > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples
  """
}




/*
 * Step 3. Run DoubletFinder
 */
process FIND_DOUBLETS {
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/finding_doublets", mode: 'copy', pattern: "${library}*"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "run_df.log" }


  container "${params.container.rsinglecell}" // todo update to one with DF

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
  tuple val(library), path("${library}_w_BCR.Rds"), emit: sobj
  path(".command.log"), emit: log

  """
  if [[ "${data_type}" == "no BCR" ]]
  then
    cp ${sobj} "${library}_w_BCR.Rds" # todo - switch to soft link
  else
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${projectDir}
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
  tuple val(library), path("${library}_w_TCR.Rds"), emit: sobj
  path(".command.log"), emit: log

  """
  if [[ "${data_type}" == "no TCR" ]]
  then
    cp ${sobj} "${library}_w_TCR.Rds" # todo - switch to soft link
  else
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${projectDir}
  fi
  
  """

}

/* 
 * Step 6. Run Seurat
 */
process SEURAT_QC {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", saveAs: { filename -> "seurat_qc.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"
  // For testing
  publishDir "${workDir}/data/single_cell_GEX/processed/${library}/automated_processing", mode: 'copy', pattern: "${library}*"

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
  path("${library}*"), emit: outfiles
  path(".command.log"), emit: log

  """
  Rscript ${projectDir}/bin/process_with_seurat_post_filter.R ${library} ${projectDir} ${raw_h5} ${params.settings.remove_demux_DBL} ${params.settings.remove_all_DBL}
  """
}

