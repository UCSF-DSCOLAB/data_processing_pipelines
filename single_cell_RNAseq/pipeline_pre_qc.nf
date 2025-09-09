nextflow.enable.dsl=2

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (workflow.success){
       println "Deleting working directory $workDir"
      // "rm -rf $workDir".execute()
    }
}


// Processes

include {
TEST_GZIP_INTEGRITY;
CELLRANGER;
CELLRANGER_VDJ;
CELLRANGER_ATAC;
FILTER_REF_VCF_ATAC;
FILTER_BARCODES;
AMULET_ATAC;
FILTER_BAM;
DSC_PILEUP;
MERGE_DSC;
FREEMUXLET_POOL;
FREEMUXLET_LIBRARY;
FMX_ASSIGN_TO_GT;
DEMUXLET_POOL;
DEMUXLET_LIBRARY;
SEPARATE_DMX;
SEPARATE_FMX;
UNMERGE_FMX;
FIND_DOUBLETS;
LOAD_SOBJ;
SEURAT_ADD_BCR;
SEURAT_ADD_TCR;
SEURAT_QC;
ARCHR_LOAD_QC
} from './modules/pipeline_tasks.nf'


// Helper functions

include {
get_c4_h5; get_c4_bam; get_c4_h5_bam_bc; 
get_c4_atac_fragments; get_c4_atac_bam; get_c4_atac_bc; get_c4_amulet_bc; get_c4_atac_peaks;
get_pool_library_meta; get_libraries_data_type_tuples;
get_pool_by_sample_count; get_library_by_sample_count; get_single_library_by_pool;
get_multi_pool_by_library ; get_library_by_pool; get_multi_library_by_pool; get_pool_vcf ; get_library_ncells;
get_vdj_tuple; get_vdj_name ; get_clonotypes; get_contigs
} from  './helpers/params_parse.nf'

include {
extractFileName
} from "./helpers/utils.nf"

log.info """\
         DSCoLab scRNASeq Pipeline
         =============================
         project_dir: ${params.project_dir}
         """
         .stripIndent()
         
workflow {
  
      // TODO: perhaps just shove everything in a bam and h5 channel, and do not differentiate between data types
      ch_gex_cite_cr = Channel.empty()
      ch_atac_cr = Channel.empty()
      ch_vdj_libs = Channel.empty()
      ch_atac_peaks = Channel.empty()

      if (params.settings.skip_cellranger){
            ch_library_info = Channel.from(get_libraries_data_type_tuples()).transpose()
              .branch{
                gex_cite: it[1] in ["GEX", "CITE"]
                bcr_tcr: it[1] in ["BCR", "TCR"]
                atac: it[1] in ["ATAC"]
            }
            ch_gex_cite_cr = ch_library_info.gex_cite // [[library, data_type, cellranger_bam, raw_h5, bc]
            .map{
              it -> [it[0], it[1], get_c4_h5(it[0]), get_c4_bam(it[0]), get_c4_h5(it[0]), get_c4_cr_filt_bc(it[0])]
            } 


            ch_vdj_libs = ch_library_info.bcr_tcr
              .map{
                it -> [it[0], it[1], get_clonotypes(it[0], it[1]), get_contigs(it[0], it[1])]
              }
              .branch { 
                        tcr: it[1].contains("TCR")
                        bcr: it[1].contains("BCR")
                    }  
            ch_atac_cr = ch_library_info.atac
              .map{
                it -> [it[0], it[1],  get_c4_atac_bam(it[0]), get_c4_atac_fragments(it[0]), get_c4_atac_bc(it[0])]
              }
            ch_atac_peaks = ch_library_info.atac.map{ it -> [it[0], get_c4_atac_peaks(it[0])]}

      } else {
            ch_library_info = Channel.from(get_libraries_data_type_tuples()).transpose()
            TEST_GZIP_INTEGRITY(ch_library_info) // -> [[library_dir, data_type]]
            
            ch_gzip_out = TEST_GZIP_INTEGRITY.out
            .branch{
              gex_cite: it[1] in ["GEX", "CITE"]
              bcr_tcr: it[1] in ["BCR", "TCR"]
              atac: it[1] in ["ATAC"]
            }

            // Run cellranger for GEX and CITE data types
            CELLRANGER(ch_gzip_out.gex_cite)
            ch_gex_cite_cr = CELLRANGER.out.bam_h5_bc // --> [[library, data_type, cellranger_bam, raw_h5, bc]]

            CELLRANGER_ATAC(ch_gzip_out.atac)
            ch_atac_cr = CELLRANGER_ATAC.out.bam_frag_bc // --> [[library, data_type, cellranger_bam, fragments, bc]]
            ch_atac_peaks = CELLRANGER_ATAC.out.peaks // --> [[library, peaks]]
            
            // Run cellranger for BCR and TCR data types
            if (params.settings.add_tcr || params.settings.add_bcr ){
                
                ch_vdj_in = ch_gzip_out.bcr_tcr.map{
                  it -> get_vdj_tuple(it[0], it[1])
                }
                CELLRANGER_VDJ(ch_vdj_in) 
                
                ch_vdj_libs = CELLRANGER_VDJ.out.vdj_csvs
                .branch{
                  tcr: it[1].contains("TCR")
                  bcr: it[1].contains("BCR")
                }
 
            }
        }
      
      // Extract all bam and h5 files
      ch_all_cr_out = ch_gex_cite_cr.concat(ch_atac_cr)
      ch_all_bam = ch_all_cr_out.map{ it -> [it[0], it[1], it[2]] } // [[library, data_type, cellranger_bam]]
      ch_all_dt_h5 = ch_all_cr_out.map{ it -> [it[0], it[1], it[3]] } // [[library, data_type, raw_h5 or frag ]]
      ch_all_bc = ch_all_cr_out.map{ it -> [it[0], it[1], it[4]] } // [[library, data_type, bc ]]
      ch_all_dt_h5_bc = ch_all_cr_out.map{ it -> [it[0], it[1], it[3], it[4]] } // [[library, data_type, raw_h5 or frag, bc ]]

      ch_atac_lib_pool = Channel.from(get_multi_library_by_pool()) // [[library, pool], [library, pool]]

      // Match libraries to pools, group by pools, and then group plp files by pools
      ch_atac_peaks_transformed = ch_atac_peaks
                                            .join(ch_atac_lib_pool)
                                            .groupTuple(by: 2)
                                            .map {it -> [it[2], it[1].flatten()]} // [pool [peak_files]]
      
      //ch_atac_lib_pool2 = Channel.from(get_multi_library_by_pool()) // [[library, pool], [library, pool]]

      FILTER_REF_VCF_ATAC(ch_atac_peaks_transformed) // [pool, ref_vcf]
      ch_atac_filt_vcf = ch_atac_lib_pool // [lib, pool]
        .map{it -> [it[1], it[0]]} // [pool, lib]
        .combine(FILTER_REF_VCF_ATAC.out.filt_vcf, by: 0) // [pool, lib, ref_vcf]
        .map{ it -> [it[1], it[2]]} // [library, ref_vcf]     
      FILTER_BARCODES(ch_all_dt_h5_bc) // [library, data_type, h5, bc ]
      ch_library_barcode = FILTER_BARCODES.out.bc_list

      //** TODO ADD IN AMULET HERE ** //
      ch_bam_bc = ch_all_bam
      .branch{
              gex_cite: it[1] in ["GEX", "CITE"]
              bcr_tcr: it[1] in ["BCR", "TCR"]
              atac: it[1] in ["ATAC"]
            }
     ch_amulet_in = ch_bam_bc.atac // [library, data_type, bam, bc]
      .join(FILTER_BARCODES.out.amulet_bc_list, by: 0)
      .map{it -> [it[0], it[2], it[3]]} // [library, bam, bc]
     ch_amulet_in.view()
     
     AMULET_ATAC(ch_amulet_in)
     ch_atac_filt_bc = AMULET_ATAC.out.filt_bc


      // Combine bam files with barcodes
      ch_library_bam_barcodes = ch_all_bam
      .join(ch_library_barcode) 
      .map{it -> [it[0], it[1], it[2], it[3], params.ref.snp_ref]} // [library, data_type, cellranger_bam, barcodes, ref_vcf]

      // Filter the bam file in prep for freemux
      FILTER_BAM(ch_library_bam_barcodes) // [library, data_type, filtered_bam]
 
      // Combine barcodes with filtered bam files
      ch_library_barcodes_filtered_bam = FILTER_BAM.out.bam_file
        .join(ch_library_barcode) // [library, data_type, filtered_bam, barcodes]
        .branch{
                gex_cite: it[1] in ["GEX", "CITE"]
                atac: it[1] in ["ATAC"]
              }
      // add snp_ref
      ch_barcodes_bam_vcf = 
        ch_library_barcodes_filtered_bam.gex_cite
          .map{it -> [it[0], it[1], it[2], it[3], params.ref.snp_ref]}
        .concat(
          ch_library_barcodes_filtered_bam.atac
          .join(ch_atac_filt_vcf)
        )

      // Run dsc_pileup
      DSC_PILEUP(ch_barcodes_bam_vcf) // [library, data_type, filtered_bam, snp_ref]
      ch_plp_files = DSC_PILEUP.out.plp_files // [library, data_type, [plp_files]]
      ch_merged_libs = Channel.empty()

     if (params.settings.merge_for_demux){
        // Only fetch multiple libraries
        ch_multi_lib_pool = Channel.from(get_multi_library_by_pool()) // [[lib, pool], [lib, pool]]
        // Match libraries to pools, group by pools, and then group plp files by pools
        ch_multi_lib_pool_transformed =  ch_plp_files.join(ch_multi_lib_pool)
                                            .groupTuple(by: 3)
                                            .map {it -> [it[3], it[1][0], it[2].flatten()]} // [pool [plp_files]]

        MERGE_DSC(ch_multi_lib_pool_transformed)
        ch_merged_libs = MERGE_DSC.out.merged_files
      }

      /*
      --------------------------------------------------------
      De-multiplex using freemuxlet or demuxlet
      --------------------------------------------------------
      */

     // We de-multiplex if we have merged libraries
     ch_sample_map = Channel.empty()
     ch_lib_vcf = Channel.empty()
     if ( params.settings.demux_method == "freemuxlet"){

        // This assumes you have at least one pool with > 1 libraries!!!
        if (params.settings.merge_for_demux) {
             // Attach the number of samples, and re-arrange input
            ch_multi_lib_transformed = ch_merged_libs
                                            .join(Channel.from(get_pool_by_sample_count()))
                                            .map{it -> [it[0], it[1], it[7], it[3], it[4], it[5], it[6]]} // [lib, data_type, num_samples, plp_files]] (excluding .tsv)
            FREEMUXLET_POOL(ch_multi_lib_transformed)

            // Combine merged files and merged tsv
            pool_tsv = ch_merged_libs.map{it -> it[0,1,2]} // [pool, data_type, pool_tsv]
            ch_freemux_transformed = pool_tsv.join(FREEMUXLET_POOL.out.merged_files)

            UNMERGE_FMX(ch_freemux_transformed)

            // TODO: I dont think this step has been fixed to handle data type!!!! ****** //
            sample_file_transformed = UNMERGE_FMX.out.samples_file
                                            .transpose() // We need to group each sub array by index [[1,2],[a,b]] -> [[1,a],[2,b]]
                                            .map {sublist ->
                                                    // Create a new sublist with the filename part and the rest of the original sublist as its own sublist
                                                    [extractFileName(sublist[0].toString()), sublist[0..-1]].flatten()
                                            }
            SEPARATE_FMX(sample_file_transformed)

            // Run freemuxlet on remaining pools with single libraries
            // Attach the number of samples, and re-arrange input
            ch_single_lib_transformed  = ch_plp_files
                                            .join(Channel.from(get_single_library_by_pool()))
                                            .join(Channel.from(get_library_by_sample_count()))
                                            .map{it -> [it[0], it[1], it[4], it[2]]} // [lib, data_type, num_of_samples, plp_files]]

            FREEMUXLET_LIBRARY(ch_single_lib_transformed)

            // appended any merged libraries
            ch_sample_map = SEPARATE_FMX.out.sample_map.mix(FREEMUXLET_LIBRARY.out.sample_map)

            ch_lib_vcf = SEPARATE_FMX.out.fmx_files.map{
              it -> [it[0], it[1], it[3]] // [library, data_type, vcf]
            }.mix(
              FREEMUXLET_LIBRARY.out.vcf
            )

        } else {
                // Run freemuxlet on all libraries, regardless if there are many libraries per pool
                // Attach the number of samples, and re-arrange input
                ch_single_lib_transformed  = ch_plp_files
                                               .join(Channel.from(get_library_by_sample_count()))
                                               .map{it -> [it[0], it[1], it[3], it[2]]} // [lib, data_type, num_of_samples, plp_files]]
                FREEMUXLET_LIBRARY(ch_single_lib_transformed)

                ch_sample_map = FREEMUXLET_LIBRARY.out.sample_map
                ch_lib_vcf =  FREEMUXLET_LIBRARY.out.vcf
          }

        } else if ( params.settings.demux_method == "demuxlet"){
             // This assumes you have at least one pool with > 1 libraries!!!
            if (params.settings.merge_for_demux) {
                ch_multi_lib_transformed = ch_merged_libs
                                                   .join(Channel.from(get_pool_vcf()))
                                                   .map{it -> [it[0], it[1], it[7], it[3], it[4], it[5], it[6]]} // [lib, data_type, vcf, plp_files]]
                DEMUXLET_POOL(ch_multi_lib_transformed) // [pool, data_type, merged_best]

                demuxlet_pool_transformed = DEMUXLET_POOL.out.merged_best
                                                  .cross(Channel.from(get_multi_pool_by_library())) // [[pool, data_type, merged_best] , [pool, lib]]
                                                  .map{it -> [it[1][1], it[0][1], it[0][2]]} // [lib, data_type, merged.best]
                
                SEPARATE_DMX(demuxlet_pool_transformed)

                // Run demuxlet on remaining single libraries
                // Attach the number of samples, and re-arrange input
                
                ch_single_lib_transformed  = ch_plp_files // [pool, data_type, files]
                                              .join(Channel.from(get_single_library_by_pool()))
                                              .map{it -> [it[0], it[3], it[1], it[2]]} // [pool, lib, data_type, files]
                                              .join(Channel.from(get_pool_vcf()))
                                              .map{it -> [it[1], it[2], it[4], it[3]]} // [lib, data_type, vcf, plp_files]]
                DEMUXLET_LIBRARY(ch_single_lib_transformed)
                // appended any merged libraries
                ch_sample_map = SEPARATE_DMX.out.sample_map.mix(DEMUXLET_LIBRARY.out.sample_map)

            } else {
                // Run demuxlet on all libraries, regardless if there are many libraries per pool
                // Attach the number of samples, and re-arrange input
   		        ch_single_lib_transformed  = Channel.from(get_pool_vcf()) // [pool, vcf]
                                                  .cross(
                                                    ch_plp_files // [library, data_type, files]
                                                    .join(Channel.from(get_library_by_pool())) // [library, data_type, files, pool]
                                                    .map{it -> [it[3], it[0], it[1], it[2]]} // [pool, lib, data_type, files]
                                                  ) // [[pool, vcf], [pool, lib, data_type, files]]
                                                  .map{it -> [it[1][1], it[0][1], it[1][2], it[1][3]]} // [lib, vcf, data_type, plp_files]]
	

                DEMUXLET_LIBRARY(ch_single_lib_transformed)
                ch_sample_map = DEMUXLET_LIBRARY.out.sample_map
            }

        }

        if (params.settings.fmx_assign_to_gt){
            ch_gt_input =  Channel.from(get_pool_vcf()) // [pool, vcf]
              .combine(Channel.from(get_library_by_pool()).map{ it -> [it[1], it[0]] }, by: 0) // [pool, vcf, lib]
              .map{it -> [it[2], it[1]]} // [lib, vcf]
              .join(ch_lib_vcf )  // [lib, ref_vcf, fmx_vcf]
              
              // TODO add checks that the reference file exists  
              FMX_ASSIGN_TO_GT(ch_gt_input) 
                    
        } 
     /*
     --------------------------------------------------------
     Set up archR
     --------------------------------------------------------
     */
     ch_dt_h5_sep = ch_all_dt_h5
            .branch{
              gex_cite: it[1] in ["GEX", "CITE"]
              bcr_tcr: it[1] in ["BCR", "TCR"]
              atac: it[1] in ["ATAC"]
            }

     ch_archr_in = ch_dt_h5_sep.atac // [library, dt, frag]
      .map{it -> [it[0], it[2]]}
     .join(ch_atac_filt_bc) 
     .join(ch_sample_map
      .map{it -> [it[0], it[2]]}
     )
     
     // [library, fragments, amulet_bc, demuxlet_out ]
     ARCHR_LOAD_QC(ch_archr_in)
}