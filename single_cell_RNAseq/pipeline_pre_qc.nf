nextflow.enable.dsl=2

// Processes

include {
TEST_GZIP_INTEGRITY;
CELLRANGER;
CELLRANGER_VDJ;
FILTER_BARCODES;
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
SEURAT_QC
} from './modules/pipeline_tasks.nf'


// Helper functions

include {
get_c4_h5; get_c4_bam; get_c4_h5_bam; get_pool_library_meta; get_libraries_data_type;
get_pool_by_sample_count; get_library_by_sample_count; get_single_library_by_pool;
get_multi_pool_by_library ; get_library_by_pool; get_multi_library_by_pool; get_pool_vcf ; get_library_ncells
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

//        /*
//         --------------------------------------------------------
//         Initial Input:
//                 - Libraries: Directory of fastqs
//         Output:
//                 - BAM / H5 files for each library
//         --------------------------------------------------------
//        */

      // TODO: perhaps just shove everything in a bam and h5 channel, and do not differentiate between data types
      ch_gex_cite_bam_h5 = Channel.empty()
      ch_bcr_tcr_bam_h5 = Channel.empty()

      if (params.settings.skip_cellranger){
            ch_gex_cite_bam_h5 =  Channel.from(get_c4_h5_bam()) // [[library, cell_ranger_bam, raw_h5]
            // TODO expand to work for VDJ as well
             ch_library_bcr_tcr = Channel.from(get_libraries_data_type()).filter { it[1] in ["BCR", "TCR"] }
                
              ch_vdj_libs = ch_library_bcr_tcr
              .map{
                it -> [it[0], it[1], get_clonotypes(it[0], it[1]), get_contigs(it[0], it[1])]
              }
              .branch { 
                        tcr: it[1].contains("TCR")
                        bcr: it[1].contains("BCR")
                    }  
            
      } else {
            library_info = get_libraries_data_type() // -> [[library_dir, data_type]]
            ch_library_info = Channel.from(library_info)
            TEST_GZIP_INTEGRITY(ch_library_info) // -> [[library_dir, data_type]]

            // Run cellranger for GEX and CITE data types
            ch_library_cite_gex = ch_library_info.filter { it[1] in ["GEX", "CITE"] }
            CELLRANGER(ch_library_cite_gex)
            ch_gex_cite_bam_h5 = CELLRANGER.out.bam_h5 // --> [[library, cell_ranger_bam, raw_h5]]

            // Run cellranger for BCR and TCR data types
            if (params.settings.add_tcr || params.settings.add_bcr ){
                ch_library_bcr_tcr = ch_library_info.filter { it[1] in ["BCR", "TCR"] }
                
                ch_vdj_in = ch_library_bcr_tcr.map{
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
    ch_all_bam = ch_gex_cite_bam_h5.mix(ch_bcr_tcr_bam_h5).map { it -> [it[0], it[1]] } // [[library, cell_ranger_bam]]
    ch_all_h5 = ch_gex_cite_bam_h5.mix(ch_bcr_tcr_bam_h5).map { it -> [it[0], it[2]] } // [[library, raw_h5 ]]

    /*
    --------------------------------------------------------
    Filter barcode list for freemuxlet
    --------------------------------------------------------
    */

    FILTER_BARCODES(ch_all_h5) // [library, barcodes]
    library_barcode = FILTER_BARCODES.out.bc_list

    /*
    --------------------------------------------------------
    Setup bam, dsc files for freemuxlet
    --------------------------------------------------------
    */

    // Combine bam files with barcodes
     ch_library_bam_barcodes = ch_all_bam.join(library_barcode) // [library, cell_ranger_bam, barcodes]

    // Filter the bam file in prep for freemux
     FILTER_BAM(ch_library_bam_barcodes) // [library, filtered_bam]

    // Combine barcodes with filtered bam files
     ch_library_barcodes_filtered_bam = library_barcode.join(FILTER_BAM.out.bam_file) // [library, barcodes, filtered_bam]

     // Run dsc_pileup
    DSC_PILEUP(ch_library_barcodes_filtered_bam) // [library, barcodes, filtered_bam]
    ch_plp_files = DSC_PILEUP.out.plp_files

     /*
     --------------------------------------------------------
     Merge multiple libraries per pool
     --------------------------------------------------------
     */

     ch_merged_libs = Channel.empty()

     if (params.settings.merge_for_demux){
        // Only fetch multiple libraries
        ch_multi_lib_pool = Channel.from(get_multi_library_by_pool()) // [[lib_dir, pool], [lib_dir, pool]]

        // Match libraries to pools, group by pools, and then group plp files by pools
        ch_multi_lib_pool_transformed = ch_plp_files
                                            .join(ch_multi_lib_pool)
                                            .groupTuple(by: 2)
                                            .map {it -> [it[2], it[1].flatten()]} // [pool [plp_files]]
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
                                            .map{it -> [it[0], it[6], it[2], it[3], it[4], it[5]]} // [lib, num_of_samples, plp_files]] (excluding .tsv)
            FREEMUXLET_POOL(ch_multi_lib_transformed)

            // Combine merged files and merged tsv
            pool_tsv = ch_merged_libs.map{it -> it[0,1]} // [pool, pool_tsv]
            ch_freemux_transformed = pool_tsv.join(FREEMUXLET_POOL.out.merged_files)

            UNMERGE_FMX(ch_freemux_transformed)
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
                                            .map{it -> [it[0], it[3], it[1]]} // [lib, num_of_samples, plp_files]]

            FREEMUXLET_LIBRARY(ch_single_lib_transformed)

            // appended any merged libraries
            ch_sample_map = SEPARATE_FMX.out.sample_map.mix(FREEMUXLET_LIBRARY.out.sample_map)

            ch_lib_vcf = SEPARATE_FMX.out.fmx_files.map{
              it -> [it[0], it[2]] // [library, vcf]
            }.mix(
              FREEMUXLET_LIBRARY.out.vcf
            )

        } else {
                // Run freemuxlet on all libraries, regardless if there are many libraries per pool
                // Attach the number of samples, and re-arrange input
                ch_single_lib_transformed  = ch_plp_files
                                               .join(Channel.from(get_library_by_sample_count()))
                                               .map{it -> [it[0], it[2], it[1]]} // [lib, num_of_samples, plp_files]]
                FREEMUXLET_LIBRARY(ch_single_lib_transformed)

                ch_sample_map = FREEMUXLET_LIBRARY.out.sample_map
                ch_lib_vcf =  FREEMUXLET_LIBRARY.out.vcf
          }

        } else if ( params.settings.demux_method == "demuxlet"){
             // This assumes you have at least one pool with > 1 libraries!!!
            if (params.settings.merge_for_demux) {
                ch_multi_lib_transformed = ch_merged_libs
                                                   .join(Channel.from(get_pool_vcf()))
                                                   .map{it -> [it[0], it[6], it[2], it[3], it[4], it[5]]} // [lib, vcf, plp_files]]
                DEMUXLET_POOL(ch_multi_lib_transformed)

                demuxlet_pool_transformed = DEMUXLET_POOL.out.merged_best
                                                  .cross(Channel.from(get_multi_pool_by_library()))
                                                  .map{it -> [it[1][1],it[0][1]]} // [lib, merged.best]

                SEPARATE_DMX(demuxlet_pool_transformed)

                // Run demuxlet on remaining single libraries
                // Attach the number of samples, and re-arrange input
                ch_single_lib_transformed  = ch_plp_files
                                              .join(Channel.from(get_single_library_by_pool()))
                                              .map{it -> [it[2], it[0], it[1]]} // [pool, lib, files]
                                              .join(Channel.from(get_pool_vcf()))
                                              .map{it -> [it[1], it[3], it[2]]} // [lib, vcf, plp_files]]
                DEMUXLET_LIBRARY(ch_single_lib_transformed)
                // appended any merged libraries
                ch_sample_map = SEPARATE_DMX.out.sample_map.mix(DEMUXLET_LIBRARY.out.sample_map)

            } else {
                // Run demuxlet on all libraries, regardless if there are many libraries per pool
                // Attach the number of samples, and re-arrange input
   		        ch_single_lib_transformed  = Channel.from(get_pool_vcf())
                                                .cross(ch_plp_files
                                                  .join(Channel.from(get_library_by_pool()))
                                                  .map{it -> [it[2], it[0], it[1]]} // [pool, lib, files]
                                                  )
                                                  .map{it -> [it[1][1], it[0][1], it[1][2]]} // [lib, vcf, plp_files]]
	

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
      Run doublet finder if specified
      --------------------------------------------------------
      */
     ch_initial_sobj = Channel.empty()
     if (params.settings.run_doubletfinder) {
        ch_doublet_input = Channel.from(get_library_ncells()).join(ch_all_h5).join(ch_sample_map) // [lib, ncells, raw_h5, fmx_cluster ]
        FIND_DOUBLETS(ch_doublet_input) // --> [lib, doublet_finder_sobj]
        ch_initial_sobj = FIND_DOUBLETS.out.sobj
     } else {
        ch_doublet_input = ch_all_h5.join(ch_sample_map) // [lib, ncells, raw_h5, fmx_cluster]
        LOAD_SOBJ(ch_doublet_input)
        ch_initial_sobj = LOAD_SOBJ.out.sobj
     }
     
    
     

     if (params.settings.add_tcr){
        SEURAT_ADD_TCR(ch_initial_sobj.join(ch_vdj_libs.tcr, by:0, remainder: true))
        ch_tcr_out = SEURAT_ADD_TCR.out.sobj
      } else {
        ch_tcr_out = ch_initial_sobj
      }

      if (params.settings.add_bcr){
        SEURAT_ADD_BCR(ch_tcr_out.join(ch_vdj_libs.bcr, by:0, remainder: true))
        ch_bcr_out = SEURAT_ADD_BCR.out.sobj
      } else {
        ch_bcr_out = ch_tcr_out
      } 

     /*
     --------------------------------------------------------
     Set up seurat object
     --------------------------------------------------------
     */
     ch_library_info = Channel.from(get_libraries_data_type()) // -> [[library_dir, data_type]]
     ch_seurat_input = ch_library_info.join(ch_bcr_out).join(ch_all_h5)
     SEURAT_QC(ch_seurat_input)

}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
