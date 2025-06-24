def date = new Date().format("yyyy-MM-dd")



process TEST_GZIP_INTEGRITY {

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "gzip_test_${date}.log" })

    input:
    tuple val(library), val(data_type)

    output:
    tuple val(library), val(data_type)

    """
    echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
    echo "[running TEST_GZIP_INTEGRITY]"

    if [[ "${data_type}"  == "ATAC" ]];
    then
      lib_to_use=${library}
      dir_use=single_nuclear_${data_type}
    else
      library_str=${library}
      dt=${data_type}
      lib_to_use=\${library_str}/"SCG"/"SC\${dt:0:1}"}
      dir_use=single_cell_${data_type}
    fi 


    echo " gzip --test ${params.project_dir}/data/${dir_use}/raw/\${lib_to_use}/\${lib_to_use}*.fastq.gz"
    gzip --test ${params.project_dir}/data/${dir_use}/raw/\${lib_to_use}/\${lib_to_use}*.fastq.gz
    
    if [[ "${data_type}" == "CITE" ]]
    then
      echo " gzip --test ${params.project_dir}/data/single_cell_GEX/raw/${library}/${library}*.fastq.gz"
      gzip --test ${params.project_dir}/data/single_cell_GEX/raw/${library}/${library}*.fastq.gz
    fi

    echo "-----------"
    """
}




/*
 * Step 1a. Run Cellranger
 */
process CELLRANGER {
  
  publishDir "${params.project_dir}/10X_library_csvs/", mode: 'copy', pattern: "${library}_libraries.csv"
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/", pattern: "cellranger/*", mode: 'copy'
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log",
   saveAs: { filename -> "cellranger_${date}.log" }


  container "${params.container.cellranger}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type)
  
  output:
  tuple val(library), val(data_type), path("cellranger/possorted_genome_bam.bam"), path("cellranger/raw_feature_bc_matrix.h5"),  path("cellranger/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: bam_h5_bc
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  """

  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running CELLRANGER]"

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
  echo " Using container ${params.container.cellranger}"
  echo " cellranger count --id=${library}  \
    --libraries=${library}_libraries.csv \
    --feature-ref=${params.ref.cite_feature_ref} \
    --transcriptome=${params.ref.transcriptome} "
  echo "-----------"

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
  publishDir "${params.project_dir}/data/single_cell_${data_type}/processed/${vdj_library}/cellranger", mode: 'copy'
  publishDir "${params.project_dir}/data/single_cell_${data_type}/logs/${vdj_library}/", mode: 'copy', 
    pattern: ".command.log", saveAs: { filename -> "cellranger__${date}.log" }

  container "${params.container.cellranger}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type), val(vdj_library) 
  
  output:
  tuple val(library), val(data_type), path("cellranger/clonotypes.csv"), path("cellranger/all_contig_annotations.csv"), emit: vdj_csvs
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  
  """

  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running CELLRANGER_VDJ]"
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
  

  echo " Using container ${params.container.cellranger}"
  echo " cellranger vdj --id="\${vdj_library}"  \
    --fastqs=\${vdj_path} \
    --reference=${params.ref.vdj_ref} \
    --chain=\${chain_type} "

  echo "-----------"


  cellranger vdj --id="\${vdj_library}"  \
    --fastqs=\${vdj_path} \
    --reference=${params.ref.vdj_ref} \
    --chain=\${chain_type} \
    --localcores=${task.cpus - 1} \
    --localmem=${task.memory.toGiga() - 2}
  
  mv ${vdj_library}/outs cellranger
  """
}

process CELLRANGER_ATAC{
  publishDir "${params.project_dir}/10X_library_csvs/", mode: 'copy', pattern: "${library}_libraries.csv"
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/", pattern: "cellranger/*", mode: 'copy'
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/", mode: 'copy', pattern: ".command.log",
   saveAs: { filename -> "cellranger_${date}.log" }


  container "${params.container.cellranger_atac}"
  containerOptions "-B ${params.ref.dir} -B ${params.project_dir} -B /scratch/"
  
  input:
  tuple val(library), val(data_type)
  
  output:
  tuple val(library), val(data_type), path("cellranger/possorted_bam.bam"), path("cellranger/fragments.tsv.gz"),  path("cellranger/singlecell.csv"), emit: bam_frag_bc
  tuple val(library), path("cellranger/peaks.bed"), emit: peaks
  path("cellranger/*"), emit: cr_out_files
  path(".command.log"), emit: log
  """

  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running CELLRANGER]"

  fastq_dir=${params.project_dir}/data/single_nuclear_ATAC/raw/${library}
  # create the config
  echo "fastqs,sample,library_type
    \${fastq_dir},${library},ATAC Seq" > ${library}_libraries.csv
  
  echo " Using container ${params.container.cellranger_atac}"
  echo "  cellranger-atac-${params.ref.cr_atac_version} count \
      --id=${library} \
      --fastqs=\${fastq_dir} \
      --sample=${library} \
      --reference=${params.ref.atac_ref} \
      --localcores=${task.cpus - 1} \
      --localmem=${task.memory.toGiga() - 2} "
  echo "-----------"

  cellranger-atac-${params.ref.cr_atac_version} count \
      --id=${library} \
      --fastqs=\${fastq_dir} \
      --sample=${library} \
      --reference=${params.ref.atac_ref} \
      --localcores=${task.cpus - 1} \
      --localmem=${task.memory.toGiga() - 2}

  mv ${library}/outs cellranger
  """
}


process FILTER_REF_VCF_ATAC{
  publishDir "${params.settings.merge_demux_dir}/${pool}", mode: 'copy', pattern: "peaks*"
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/logs/${pool}/", mode: 'copy', pattern: ".command.log",
   saveAs: { filename -> "filter_atac_vcf_${date}.log" }

  container "${params.container.bed_bcftools}"
  containerOptions "-B ${params.ref.dir}" 

  input:
    tuple val(pool), val(peak_files)
  output:
    tuple val(pool), path("peaks.vcf.gz"), emit: filt_vcf
    tuple val(pool), path("peaks.all.bed"), emit: merged_peaks
    path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FILTER_REF_VCF_ATAC]"
  echo " using container ${params.container.bed_bcftools}"

  touch peaks.all.bed

  for f in "${peak_files[@]}"; do
      bedtools sort -i ${f} | \
      bedtools merge -i - | cut -f 1,2,3 >> peaks.all.bed
  done

  bcftools view -R peaks.all.bed -O z -o peaks.vcf.gz ${params.ref.snp_ref}

  """
}


/*
 * Step 2a. Create barcode list
 */ 
process FILTER_BARCODES{
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet", mode: 'copy', pattern: "barcodes_of_interest.filt.list"
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", saveAs: { filename -> "filter_bc_${date}.log" }
  
  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: "barcodes_of_interest.filt.list" )

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "filter_bc_${date}.log" })

  container "${params.container.rsinglecell}"
  
  input:
  tuple val(library), val(data_type), path(raw_h5), path(bc) 
  
  output:
  tuple val(library), path("barcodes_of_interest.filt.list"), emit: bc_list
  path(".command.log"), emit: log
  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FILTER_BARCODES]"
  echo " using container ${params.container.rsinglecell}"

  if [ "${data_type}"== "ATAC" ];
  then
    echo " Rscript ${projectDir}/bin/make_valid_barcodelist_atac.R ${bc} ${params.settings.minfeature}"  
    echo "-----------"
    Rscript ${projectDir}/bin/make_valid_barcodelist_atac.R ${bc} ${params.settings.minfeature}
  else
    echo " Rscript ${projectDir}/bin/make_valid_barcodelist.R ${raw_h5} ${params.settings.minfeature} ${params.settings.mincell}"  
    echo "-----------"
    Rscript ${projectDir}/bin/make_valid_barcodelist.R ${raw_h5} ${params.settings.minfeature} ${params.settings.mincell}
  fi

  
  """
}

process AMULET_ATAC {
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/", mode: 'copy', pattern: ".command.log",
   saveAs: { filename -> "amulet_${date}.log" }
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/", mode: 'copy', pattern: "amulet/*"
 
  container "${params.container.amulet}"
  containerOptions "-B ${params.ref.dir}" 

  input:
  tuple val(library), path(bc), path(bam)

  output:
  tuple val(library), path("post_amulet_barcodes_of_interest.filt.list"), emit: filt_bc
  path(".command.log"), emit: log

  """

  ${SOFTWARE_DIR}/AMULET/AMULET.sh --forcesorted --bcidx 0 --cellidx 0 --iscellidx 9 ${bam} ${bc} ${params.ref.chr_list} \
    ${params.ref.atac_blacklist} amulet ${SOFTWARE_DIR}/AMULET/

  Rscript filter_barcodes_amulet.R ${bc} 

  """
}




/*
 * alternately, filter barcodes based on QC cutoffs
 */
process SEURAT_PRE_FMX_QC {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> "pre_fmx_qc_${date}.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter", mode: 'copy', pattern: "${library}*"

  container "${params.container.rsinglecell}" 
  containerOptions "-B ${params.settings.default_qc_cuts_dir}"
  
  input:
  tuple val(library), val(data_type), path(raw_h5), path(cr_filt_bc)

  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library), path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}*"), emit: outfiles
  path("${library}_quantiles_pre.tsv"), emit: quantiles
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_PRE_FMX_QC]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} "null" ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}"
  echo " Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${data_type} ${library}_initial_raw.rds ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}"
  echo "-----------"

  Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} "null" ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}
  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${data_type} ${library}_initial_raw.rds ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}
  
  """
}




process SEURAT_PRE_FMX_FILTER {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> "pre_fmx_filter_${date}.log" }
  publishDir "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter", mode: 'copy', pattern: "${library}*" 

  publishDir(
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing" }, 
            mode: 'copy', pattern: "${library}_cutoffs.csv")
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
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_PRE_FILTER]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/filter_barcodes_pre_fmx.R ${library} ${projectDir}"
  echo "-----------"

  Rscript ${projectDir}/bin/filter_barcodes_pre_fmx.R ${library} ${projectDir} 
  
  """
}



/*
 * Step 2b. Filter bam file
 */
process FILTER_BAM {
  container "${params.container.popscle}"
  containerOptions "-B ${params.ref.fmx_dir}"

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "filter_bam_${date}.log" })

  input:
  tuple val(library), val(data_type), path(cr_bam), path(barcodes), path(ref_vcf) 

  output:
  tuple val(library), val(data_type), path("${library}_filtered.bam"), emit: bam_file
  path(".command.log"), emit: log
  
  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FILTER_BAM]"
  echo " using container ${params.container.popscle}"
  echo " bash filter_bam_for_dsc_pileup.sh \
    ${cr_bam} \
    ${barcodes} \
    ${ref_vcf} \
    ${library}_filtered.bam"
  echo "-----------"

  bash filter_bam_for_dsc_pileup.sh \
    ${cr_bam} \
    ${barcodes} \
    ${ref_vcf} \
    ${library}_filtered.bam
  
  """
}

/*
 * Step 2c. Run DSC-pileup
 */
process DSC_PILEUP{
  container "${params.container.popscle}"
  containerOptions "-B ${params.ref.fmx_dir}"

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "dsc_pileup_${date}.log" })
  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/freemuxlet" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet" }, 
            mode: 'copy', pattern: "${library}*.gz")

  input:
  tuple val(library), val(data_type), path(filtered_bam), path(barcodes), path(snp_ref)

  output:
  tuple val(library), val(data_type), path("${library}*.gz"), emit: plp_files
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running DSC-PILEUP]"
  echo " using container ${params.container.popscle}"
  echo " popscle dsc-pileup --sam ${filtered_bam} \
                   --tag-group CB \
                   --tag-UMI UB \
                   --vcf ${snp_ref} \
                   --group-list ${barcodes} \
                   --out ${library}"
  echo "-----------"


  popscle dsc-pileup --sam ${filtered_bam} \
                   --tag-group CB \
                   --tag-UMI UB \
                   --vcf ${snp_ref} \
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

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${pool}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${pool}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "merged_fmx_${date}.log" })

  container "${params.container.python}"

  input:
  // TODO: expand out what these path files are...
  tuple val(pool), val(data_type), path(pool_files)
  
  output:
  tuple val(pool), val(data_type), path("${pool}.tsv"), path("merged.plp.gz"), path("merged.var.gz"), path("merged.cel.gz"), path("${pool}.barcodes.gz"), emit: merged_files
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running MERGE_DSC]"
  echo " using container ${params.container.python}"
  echo " python ${projectDir}/bin/merge_freemuxlet_dsc_pileups.py --freemuxlet_dir_tsv ${pool}.tsv --ignore_diff_lengths_error"
  echo "-----------"

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

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${pool}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${pool}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "freemuxlet_${date}.log" })

  container "${params.container.popscle}"

  input:
  tuple val(pool), val(data_type), val(nsamples), path(merged_plp), path(merged_var), path(merged_cel), path(merged_barcodes)
  
  output:
  tuple val(pool), val(data_type), path("merged.clust1.samples.gz"), path("merged.lmix"), path('merged.clust1.vcf.gz'), emit: merged_files
  tuple val(pool), val(data_type), path("merged.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FREEMUXLET_POOL]"
  echo " using container ${params.container.popscle}"
  echo " popscle freemuxlet --plp merged \
                     --out merged \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed} \
                     --group-list ${merged_barcodes}"
  echo "-----------"


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
  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "freemuxlet_${date}.log" })

  publishDir(
    path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/" }, 
            mode: 'copy', pattern: "${library}*")

  container "${params.container.popscle}"

  input:
  tuple val(library), val(data_type), val(nsamples), path(plp_files)
  
  output:
  tuple val(library), val(data_type), path("${library}*"), emit: fmx_files
  tuple val(library), val(data_type), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
  tuple val(library), val(data_type), path("${library}.clust1.vcf.gz"), emit: vcf
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FREEMUXLET_LIBRARY]"
  echo " using container ${params.container.popscle}"
  echo " popscle freemuxlet --plp ${library} \
                     --out ${library} \
                     --nsample ${nsamples} \
                     --seed ${params.settings.randomseed}"
  echo "-----------"


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
 * Run assign to gt for a pool or library 
 */
process FMX_ASSIGN_TO_GT {
  publishDir "${params.project_dir}/fmx_assign_to_gt/${pool}/", mode: 'copy', pattern: "${pool}_gtcheck*"
  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "fmx_assign_to_gt_${date}.log" })

  container "${params.container.rplus_bcftools}"
  containerOptions "-B ${params.project_dir} -B ${params.settings.ref_vcf_dir}"

  input:
  tuple val(pool), val(data_type), val(ref_vcf), path(fmx_vcf) 

  output:
  tuple val(pool), val(data_type), path("${pool}*"), emit: outfiles
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
  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/demuxlet" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/demuxlet" }, 
            mode: 'copy', pattern: "${library}*")

  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "demuxlet_${date}.log" })
  container "${params.container.popscle}"
  
  input:
  tuple val(library), val(data_type), path(vcf), path(plp_files)

  output:
  tuple val(library), val(data_type), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
  path("${library}*"), emit: dmx_files
  path(".command.log"), emit: log
  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running DEMUXLET_LIBRARY]"
  echo " using container ${params.container.popscle}"
  echo " popscle demuxlet --plp ${library} \
                     --out ${library} \
                     --vcf ${vcf} \
                     --field GT"
  echo "-----------"
  
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

  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "demuxlet_${date}.log" })
  container "${params.container.popscle}"
  
  input:
  tuple val(pool), val(data_type), path(vcf), path(merged_plp), path(merged_var), path(merged_cel), path(merged_barcodes)

  output:
  tuple val(pool), val(data_type), path("merged.best"), emit: merged_best
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running DEMUXLET_POOL]"
  echo " using container ${params.container.popscle}"
  echo " popscle demuxlet --plp merged \
                     --out merged \
                     --vcf ${vcf} \
                     --field GT"
  echo "-----------"

  popscle demuxlet --plp merged \
                     --out merged \
                     --vcf ${vcf} \
                     --field GT 
  
   """
}

process SEPARATE_DMX {

  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/demuxlet" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/demuxlet" }, 
            mode: 'copy', pattern: "${library}*")

  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "separate_dmx_${date}.log" })


  input: 
   tuple val(library), val(data_type), path(merged_best)
  output:
   tuple val(library), val(data_type), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path("${library}.clust1.samples.gz"), emit: full_sample_file
   path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEPARATE_DMX]"
  echo " library: ${library}"
  echo " merged_best: ${merged_best}"
  echo "-----------"

  head -1 ${merged_best} > ${library}.clust1.samples0
  grep "${library}\\s" ${merged_best} >> ${library}.clust1.samples0
  sed \"s/--${library}//g\" ${library}.clust1.samples0 > ${library}.clust1.samples
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} "${library}.clust1.samples" > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples
  """
}





process UNMERGE_FMX {
  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "unmerge_fmx_${date}.log" })

  container "${params.container.python}"

  input:
  tuple val(pool), val(data_type), path(fmx_tsv), path(merged_samples), path(merged_lmix), path(merged_vcf)
  
  output:
  tuple path("*vcf.gz"), path("*samples.gz"), path("*.lmix"), emit: samples_file
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running UNMERGE_FMX]"
  echo " using container ${params.container.python}"
  echo " python ${projectDir}/bin/unmerge_freemuxlet_dsc_pileups.py --freemuxlet_dir_tsv ${fmx_tsv}"
  echo "-----------"

  # unmerge
  python ${projectDir}/bin/unmerge_freemuxlet_dsc_pileups.py --freemuxlet_dir_tsv ${fmx_tsv} 
  rm merged*
  """
}

process SEPARATE_FMX {
  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/freemuxlet" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet" }, 
            mode: 'copy', pattern: "${library}*")  
  publishDir ( path: { "${data_type}" == "ATAC" ? 
            "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/" : 
            "${params.project_dir}/data/single_cell_GEX/logs/${library}/" }, 
            mode: 'copy', pattern: ".command.log",
       saveAs: { filename -> "separate_fmx_${date}.log" })
  input:
   tuple val(library), val(data_type), path(vcf_file), path(sample_file), path(lmix_file)
  output:
   tuple val(library), val(data_type), path("${library}.clust1.samples.gz"), path("${library}.clust1.vcf.gz"), path("${library}.lmix"), emit: fmx_files
   tuple val(library), val(data_type), path("${library}.clust1.samples.reduced.tsv"), emit: sample_map
   path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEPARATE_FMX]"
  echo " library: ${library}"
  echo "-----------"

  gunzip -f ${library}.clust1.samples.gz
  awk {'printf (\"%s\t%s\t%s\t%s\t%s\\n\", \$2, \$3, \$4, \$5, \$6)'} ${library}.clust1.samples > ${library}.clust1.samples.reduced.tsv
  gzip -f ${library}.clust1.samples
  """
}


/*
 * Step 3. Run DoubletFinder
 */
process FIND_DOUBLETS {
  publishDir ( 
    path: "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> 
          params.settings.demux_method.equals("demuxlet") ? "run_df_dmx_${date}.log" : "run_df_${date}.log" }
  )
  publishDir( 
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/finding_doublets_dmx" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/finding_doublets" }, 
            mode: 'copy', pattern: "${library}*"
  )
  container "${params.container.rsinglecell}" 

  input:
  tuple val(library), val(ncells_loaded), path(raw_h5), path(fmx_clusters)

  output:
  tuple val(library), path("${library}_seurat_object_findingDoublets.rds"), emit: sobj
  path("${library}*"), emit: df_stats
  path(".command.log"), emit: log

  """
  arr="${params.settings}"
  if [[ \${arr} == *"use_inter_dbl_rate"* ]];
  then
    use_inter_dbl_rate=${params.settings.use_inter_dbl_rate}
  else
    echo "'use_inter_dbl_rate' parameter not found, defaulting to true"
    use_inter_dbl_rate="true"
  fi

  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running FIND_DOUBLETS]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/find_doublets.R ${raw_h5} ${fmx_clusters} ${library} ${ncells_loaded} ${params.settings.minfeature} ${params.settings.mincell} ${params.settings.randomseed} \${use_inter_dbl_rate} ${projectDir}"
  echo "-----------"
  
  Rscript ${projectDir}/bin/find_doublets.R ${raw_h5} ${fmx_clusters} ${library} ${ncells_loaded} ${params.settings.minfeature} ${params.settings.mincell} ${params.settings.randomseed} \${use_inter_dbl_rate} ${projectDir}
  
  """
}

process LOAD_SOBJ {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", 
    saveAs: { filename -> "load_sobj_${date}.log" }

  container "${params.container.rsinglecell}" // todo update to one with DF

  input:
  tuple val(library), path(raw_h5), path(fmx_clusters)

  output:
  tuple val(library), path("${library}*.rds"), emit: sobj
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running LOAD_SOBJ]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} ${fmx_clusters} ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}"
  echo "-----------"
  
  Rscript ${projectDir}/bin/load_sobj.R ${raw_h5} ${fmx_clusters} ${library} ${params.settings.minfeature} ${params.settings.mincell} ${projectDir}
  
  """
}








/* 
 * Step 5a. Seurat add vdj
 */
process SEURAT_ADD_BCR {
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", 
    saveAs: { filename -> "add_bcr_${date}.log" }

  container "${params.container.rsinglecell}"

  
  input:
  tuple val(library), path(sobj), val(data_type), path(clonotypes_csv), path(contig_csv)
  
  output:
  tuple val(library), path("${library}_w_BCR.RDS"), emit: sobj
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_ADD_BCR]"
  echo " using container ${params.container.rsinglecell}"
  echo " data_type: ${data_type}"

  if [[ "${data_type}" == "no BCR" ]]
  then
    cp ${sobj} "${library}_w_BCR.RDS" # todo - switch to soft link
  else
    echo " Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${clonotypes_csv} ${contig_csv} ${projectDir}
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${clonotypes_csv} ${contig_csv} ${projectDir}
  fi
  
  echo "-----------"
  """
}

/* 
 * Step 5b. Seurat add vdj
 */
process SEURAT_ADD_TCR {
    
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", 
    saveAs: { filename -> "add_tcr_${date}.log" }

  container "${params.container.rsinglecell}"

  
  input:
  tuple val(library), path(sobj), val(data_type), path(clonotypes_csv), path(contig_csv)
  
  output:
  tuple val(library), path("${library}_w_TCR.RDS"), emit: sobj
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_ADD_TCR]"
  echo " using container ${params.container.rsinglecell}"
  echo " data_type: ${data_type}"

  if [[ "${data_type}" == "no TCR" ]]
  then
    cp ${sobj} "${library}_w_TCR.RDS" # todo - switch to soft link
  else
    echo "Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${clonotypes_csv} ${contig_csv} ${projectDir}"
    Rscript ${projectDir}/bin/seurat_add_vdj.R ${library} ${sobj} ${data_type} ${clonotypes_csv} ${contig_csv} ${projectDir}
  fi
    
  echo "-----------"
  """
}

/* 
 * Step 6. Run Seurat
 */
process SEURAT_QC {
  publishDir ("${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> params.settings.demux_method.equals("demuxlet") ? "seurat_qc_dmx_${date}.log" : "seurat_qc_${date}.log" })
  publishDir( 
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing" }, 
            mode: 'copy', pattern: "${library}*"
  )

  container "${params.container.rsinglecell}"
  containerOptions "-B ${params.settings.default_qc_cuts_dir}"

  input:
  tuple val(library), val(main_dt), path(doublet_finder_sobj), path(raw_h5), path(cr_filt_bc)
  
  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}*"), emit: outfiles
  path("${library}_quantiles_pre.tsv"), emit: quantiles
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_QC]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}"
  echo "-----------"

  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}
  """
}

/* 
 * Step 6b. Run Seurat
 */
process SEURAT_LOAD_POST_QC {
  publishDir ("${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> params.settings.demux_method.equals("demuxlet") ? "seurat_qc_dmx_${date}.log" : "seurat_qc_${date}.log" })
  publishDir( 
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing" }, 
            mode: 'copy', pattern: "${library}*"
  )

  // For testing
  publishDir( 
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${workDir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx" : 
            "${workDir}/data/single_cell_GEX/processed/${library}/automated_processing" }, 
            mode: 'copy', pattern: "${library}*"
  )

  container "${params.container.rsinglecell}"

  input:
  tuple val(library), val(main_dt), path(doublet_finder_sobj), path(raw_h5), path(cr_filt_bc), path(cutoffs)
  
  output:
  tuple val(library), path("${library}_raw.rds"), emit: qc_output
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path(".command.log"), emit: log

  """
  echo "Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}"
  Rscript ${projectDir}/bin/process_with_seurat.R ${library} ${main_dt} ${doublet_finder_sobj} ${projectDir} ${params.settings.default_qc_cuts_dir}/${params.settings.default_qc_cuts_file} ${raw_h5} ${cr_filt_bc}
  """
}


/* 
 * Step 5b. Run Post filter
 */
process SEURAT_POST_FILTER {
  publishDir (
    path: { params.settings.demux_method.equals("demuxlet") ? 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx" : 
            "${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing" },  
            mode: 'copy', pattern: "${library}*"
  )
  publishDir "${params.project_dir}/data/single_cell_GEX/logs/${library}/", mode: 'copy', pattern: ".command.log", 
    saveAs: { filename -> params.settings.demux_method.equals("demuxlet") ? "post_filter_dmx_${date}.log" : "post_filter_${date}.log" }
    
  container "${params.container.rsinglecell}"

  input:
  tuple val(library), path(cutoffs), path(raw_sobj), path(raw_h5)
  
  output:
  tuple val(library), path("${library}_filtered.rds"), emit: post_process
  tuple val(library),path("${library}_cutoffs.csv"), emit: cutoffs_file
  path("${library}*"), emit: outfiles
  path(".command.log"), emit: log

  """
  echo "[\$(date '+%d/%m/%Y %H:%M:%S')]"
  echo "[running SEURAT_POST_FILTER]"
  echo " using container ${params.container.rsinglecell}"
  echo " Rscript ${projectDir}/bin/process_with_seurat_post_filter.R ${library} ${projectDir} ${raw_h5} ${params.settings.remove_demux_DBL} ${params.settings.remove_all_DBL}"
  echo "-----------"

  Rscript ${projectDir}/bin/process_with_seurat_post_filter.R ${library} ${projectDir} ${raw_h5} ${params.settings.remove_demux_DBL} ${params.settings.remove_all_DBL}
  """
}

process ARCHR_LOAD_QC {
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/", pattern: "archR/*", mode: 'copy'

  publishDir "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", 
    saveAs: { filename -> "archR_load_qc${date}.log" }

  container "${params.container.rsinglecell}"


  input:
  tuple val(library), path(fragments), path(amulet_bc), path(sample_map)

  output:
  tuple val(library), path(project_rdata), emit: rdata
  path("archR*"), emit: outfiles
  path(".command.log"), emit: log

  """

  Rscript archR_load_qc.R ${library} ${fragments} ${amulet_bc} ${sample_map}

  """
}

process ARCHR_POST_QC {
  publishDir "${params.project_dir}/data/single_nuclear_ATAC/processed/${library}/", pattern: "archR/*", mode: 'copy'

  publishDir "${params.project_dir}/data/single_nuclear_ATAC/logs/${library}/", mode: 'copy', 
    pattern: ".command.log", 
    saveAs: { filename -> "archR_post_qc_${date}.log" }

  container "${params.container.rsinglecell}"


  input:
  tuple val(library), path(project_rdata)

  output:
  tuple val(library), path(project_rdata), path("${library}_cutoffs.csv"), emit: rdata_cuts
  path("archR*"), emit: outfiles

  path(".command.log"), emit: log

  """

  Rscript archR_post_qc.R ${library} ${project_rdata}

  """
}
