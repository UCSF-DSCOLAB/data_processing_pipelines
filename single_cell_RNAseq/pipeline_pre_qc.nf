nextflow.enable.dsl=2

// This pipeline expects the following to be provided in json form on the command line
/* params.project_dir  
 * params.pools is a list of pools
 *  - pool:
 *      - nsamples
 *      - vcf  // optional
 *      - libraries:
 *        - library1: ncells_loaded
 *        - library2: ncells_loaded
 *      - data_types: [ ] // anyOf: CITE, TCR, BCR
 * see "nextflow_schema.json" for more details
 */

log.info """\
         DSCoLab scRNASeq Pipeline
         =============================
         project_dir: ${params.project_dir}
         pools : ${params.pools}
         """
         .stripIndent()


// Param Parsing

def get_c4_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}
def get_c4_bam(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/possorted_genome_bam.bam", checkIfExists: true)
}

def get_c4_h5_bam(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.dir, get_c4_bam(library.dir), get_c4_h5(library.dir)]
           }
        }
}

def get_libraries_data_type(){
    return params.pools.collectMany {
                pool -> pool.libraries.collect {
                    library -> [library.dir, library.data_types.join(",")]
                }
           }
}


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

    ch_all_bam = ch_gex_cite_bam_h5.mix(ch_bcr_tcr_bam_h5).map { it -> [it[0], it[1]] } // [[library, cell_ranger_bam ]]
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
     DSC_PILEUP(ch_library_barcodes_filtered_bam) // [library, barcodes, filtered_bam] )
     ch_plp_files = DSC_PILEUP.out.plp_files

 }


//      /*
//      --------------------------------------------------------
//      Run freemuxlet
//      --------------------------------------------------------
//      */
//
//      ch_sample_map = Channel.empty()
//
//      // ---- merge freemuxlet by pool if specified ---- //
//      if (params.settings.merge_for_demux){
//
//         // prepare input channel: first group by pool
//         ch_dsc_by_pool = ch_plp_files
//           .combine(ch_lib_pools, by:0) // add pool
//           .map { it ->  it[2,0,1] } // pool, lib, list_f
//           .groupTuple() // group by pool
//           .combine(ch_pool_nlibs, by:0) // add nlibs
//           .combine(ch_pool_ns, by:0) // add num samples for freemux
//           .combine(ch_pool_vcf, by:0) // add vcf for demux --> [pool, libs, list_f, nlibs, nsamples, vcf]
//           .map {it -> it[0,3,1,4,5,2]} // --> [pool, nlibs, libs, nsamples, vcf, list_f]
//
//         // second separate out the ones to merge (multi-lib pool) & not (single lib pool)
//         ch_merge_in = ch_dsc_by_pool
//           .branch{
//             merge: it[1] > 1 // multi-lib
//               return tuple(it[0], it[3], it[4], it[5].flatten().collect()) // pool, nsamples, vcf, files
//             skip_merge: it[1] == 1 // single lib
//               return tuple(it[2][0], it[3], it[4], it[5][0]) // libs, nsamples, vcf, files
//          }
//
//
//         // run merge on the multi-lib ones
//         // first separate for dsc vs fmx steps
//         ch_merge_in_split = ch_merge_in.merge.multiMap{ it ->
//           dsc_merge_in: it[0,3] // pool, files
//           pool_ns_vcf: it[0,1,2] // pool, nsamples, vcf
//         }
//
//         MERGE_DSC(ch_merge_in_split.dsc_merge_in) // --> [pool1, pool1.tsv, merged_plp, merged_cl, merged_var, barcodes]
//
//         ch_merged = MERGE_DSC.out.merged_files.multiMap{ it ->
//           demux_in: it[0,2,3,4,5] // [pool1, merged_plp, merged_cl, merged_var, barcodes]
//           unmerge_in: it[0,1]  // [pool, pool_tsv]
//         }
//
//         if (params.settings.demux_method == "freemuxlet"){
//           // run single library FMX on the single-lib ones
//           FREEMUXLET_LIBRARY(ch_merge_in.skip_merge.map{
//             it -> it[0,1,3] // remove VCF
//             })
//
//           // run freemuxlet on merged
//           ch_fmx_in = ch_merge_in_split.pool_ns_vcf
//             .map{ it -> it[0,1]} // remove vcf
//             .combine(ch_merged.demux_in, by: 0) // --> [pool1,  numsamples, ...]
//           FREEMUXLET_POOL(ch_fmx_in)
//
//           // separate demulitplexing output files
//           ch_unmerge_in = ch_merged.unmerge_in.combine( FREEMUXLET_POOL.out.merged_files, by: 0)
//           UNMERGE_FMX(ch_unmerge_in) // --> [pool, list_of_f]
//
//           ch_separate = UNMERGE_FMX.out.samples_file
//             .flatten()
//             .map { it ->
//             def lib = it.name.toString().tokenize('.').get(0) // extracts the library
//             return tuple(lib, it)
//             }
//             .groupTuple() // group by library
//
//           SEPARATE_FMX(ch_separate) // --> lib, fmx_clus
//           ch_sample_map = SEPARATE_FMX.out.sample_map
//             .mix(FREEMUXLET_LIBRARY.out.sample_map)
//         }
//         if (params.settings.demux_method == "demuxlet") {
//           DEMUXLET_LIBRARY(ch_merge_in.skip_merge.map{
//             it -> it[0,2,3] // remove nsamples
//             })
//
//           // run freemuxlet on merged
//           ch_dmx_in = ch_merge_in_split.pool_ns_vcf
//             .map{ it -> it[0,2]} // remove nsamples
//             .combine(ch_merged.demux_in, by: 0) // --> [pool1,  vcf, ...]
//           DEMUXLET_POOL(ch_dmx_in)
//
//           ch_separate = DEMUXLET_POOL.out.merged_best
//             .combine(ch_pool_libs, by: 0)
//             .map{ it[2,1] } // get rid of pool
//
//           SEPARATE_DMX(ch_separate)
//           ch_sample_map = SEPARATE_DMX.out.sample_map
//             .mix(DEMUXLET_LIBRARY.out.sample_map)
//         }
//       }
//
//        // ---- if not merging, run on each library individually --- //
//        else {
//          if (params.settings.demux_method == "freemuxlet"){
//           ch_fmx_in = ch_plp_files
//               .combine(ch_lib_pools, by: 0)
//               .map{it -> it[2,0,1]}
//               .combine(ch_pool_ns, by:0) // add numsamples
//               .map {it -> it[1,3,2]} //  --> [library, nsamples, files]
//           FREEMUXLET_LIBRARY(ch_fmx_in)
//           ch_sample_map = FREEMUXLET_LIBRARY.out.sample_map
//          }
//          if (params.settings.demux_method == "demuxlet"){
//           ch_dmx_in = ch_plp_files
//               .combine(ch_lib_pools, by: 0)
//               .map{it -> it[2,0,1]}
//               .combine(ch_pool_vcf, by:0) // add vcf
//               .map {it -> it[1,3,2]} //  --> [library, vcf, files]
//
//           DEMUXLET_LIBRARY(ch_dmx_in)
//           ch_sample_map = DEMUXLET_LIBRARY.out.sample_map
//          }
//
//       }
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
