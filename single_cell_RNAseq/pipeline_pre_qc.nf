nextflow.enable.dsl=2

log.info """\
         DSCoLab scRNASeq Pipeline
         =============================
         project_dir: ${params.project_dir}
         pools : ${params.pools}
         """
         .stripIndent()


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
DEMUXLET_POOL;
DEMUXLET_LIBRARY;
SEPARATE_DMX;
UNMERGE_FMX;
SEPARATE_FMX;
FIND_DOUBLETS;
LOAD_SOBJ;
SEURAT_ADD_BCR;
SEURAT_ADD_TCR;
SEURAT_QC
} from './modules/pipeline_tasks.nf'


// Helper functions

include {
get_c4_h5; get_c4_bam; get_c4_h5_bam; get_pool_library_meta; get_libraries_data_type; get_pools_with_multi_library;
get_pool_by_sample_count; get_library_by_sample_count; get_single_library_by_pool; get_multi_library_by_pool ; get_pool_vcf
} from  './helpers/params_parse.nf'

include {
extractFileName
} from "./helpers/utils.nf"

workflow {

       /*
        --------------------------------------------------------
        STEP 1
              Input:
                - Libraries: Directory of fastqs
              Output:
                - BAM / H5 files for each library
        --------------------------------------------------------
       */

      // TODO: perhaps just shove everything in a bam and h5 channel, and do not differentiate between data types
      ch_gex_cite_bam_h5 = Channel.empty()
      ch_bcr_tcr_bam_h5 = Channel.empty()

      if (params.settings.skip_cellranger){
            ch_gex_cite_bam_h5 =  Channel.from(get_c4_h5_bam()) // [[library, cell_ranger_bam, raw_h5]
            // TODO expand to work for VDJ as well
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
                CELLRANGER_VDJ(ch_library_bcr_tcr)
                ch_bcr_tcr_bam_h5 = CELLRANGER_VDJ.out.bam_h5
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

     if ( params.settings.demux_method == "freemuxlet"){

        // Did we merge any libraries?
        ch_sample_map_merged = Channel.empty()
        if (params.settings.merge_for_demux) {

             // Run freemuxlet on merged libraries
             // Attach the number of samples, and re-arrange input
            ch_multi_lib_transformed = ch_merged_libs
                                            .join(Channel.from(get_pools_with_multi_library()))
                                            .join(Channel.from(get_pool_by_sample_count()))
                                            .map{it -> [it[0], it[6], it[2], it[3], it[4], it[5]]} // [lib, num_of_samples, plp_files]] (excluding .tsv)
            FREEMUXLET_POOL(ch_multi_lib_transformed)

            // Combine merged files and merged tsv
            pool_tsv = ch_merged_libs.map{it -> it[0,1]} // [pool, pool_tsv]
            ch_freemux_transformed = pool_tsv.join(FREEMUXLET_POOL.out.merged_files)

            UNMERGE_FMX(ch_freemux_transformed)
            sample_file_transformed = UNMERGE_FMX.out.samples_file
                                            .transpose() // We need to group each sub array by index [[1,2],[a,b]] -> [[1,a],[2,b]]
                                            .collect {
                                                sublist ->
                                                    // Create a new sublist with the filename part and the rest of the original sublist as its own sublistu
                                                    return [
                                                        [extractFileName(sublist[0].toString()), sublist[0..-1]]
                                                    ]
                                            }
            sample_file_transformed.view()
            SEPARATE_FMX(sample_file_transformed.filter{it -> it[0] =='TEST-POOL-DM1-SCG1'})
            //ch_sample_map_merged = SEPARATE_FMX.out.sample_map


        }

        // Run freemuxlet on single libraries
        // Attach the number of samples, and re-arrange input
         ch_single_lib_transformed  = ch_plp_files
                                               .join(Channel.from(get_single_library_by_pool()))
                                               .join(Channel.from(get_library_by_sample_count()))
                                               .map{it -> [it[0], it[3], it[1]]} // [lib, num_of_samples, plp_files]]

        FREEMUXLET_LIBRARY(ch_single_lib_transformed)

        // appended any merged libraries
        ch_sample_map = ch_sample_map_merged.mix(FREEMUXLET_LIBRARY.out.sample_map)

        } else if ( params.settings.demux_method == "demuxlet"){


            if (!params.settings.merge_for_demux) {


                }

        }





//
//       /*
//       --------------------------------------------------------
//       Run doublet finder if specified
//       --------------------------------------------------------
//       */
//       ch_initial_sobj = Channel.empty()
//       ch_df_in = ch_cr_out.df_in.combine(ch_sample_map, by: 0) // --> [lib, ncells, raw_h5, fmx_clusters]
//       if (params.settings.run_doubletfinder) {
//         ch_df_in2 = ch_lib_ncells.combine(ch_df_in, by: 0)
//         FIND_DOUBLETS(ch_df_in2) // --> [lib, doublet_finder_sobj]
//         ch_initial_sobj = FIND_DOUBLETS.out.sobj
//       }
//       else {
//         LOAD_SOBJ(ch_df_in)
//         ch_initial_sobj = LOAD_SOBJ.out.sobj
//       }
//
//       /*
//       --------------------------------------------------------
//       Set up seurat object
//       --------------------------------------------------------
//       */
//       // add TCR & BCR data
//       ch_tcr_out = Channel.empty()
//       ch_bcr_out = Channel.empty()
//       ch_tcr_out = (params.settings.add_tcr) ? SEURAT_ADD_TCR(ch_initial_sobj.combine(ch_vdj_libs.tcr, by:0)).out.sobj : ch_initial_sobj
//       ch_bcr_out = (params.settings.add_bcr) ? SEURAT_ADD_BCR(ch_tcr_out.combine(ch_vdj_libs.bcr, by:0)).out.sobj : ch_tcr_out
//
//       // set up seurat QC
//       ch_seurat_qc_in = ch_library_dt.seurat_in
//       .combine(ch_bcr_out, by:0)
//       .combine(ch_cr_out.seurat_in, by:0) // <-- [lib, data_type, doublet_finder_sobj, raw_h5]
//
//       SEURAT_QC(ch_seurat_qc_in)
//
// }
//
// workflow.onComplete {
//     log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
// }
//
//
//
//
//
//
//
}