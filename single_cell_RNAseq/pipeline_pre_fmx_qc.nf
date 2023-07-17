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



/* 
 * define getters for required tasks
 */

def get_libraries (pool) {
  return params.pools[pool].libraries.keySet()
}
def get_data_types (pool, library){
  return params.pools[pool].libraries[library].data_types
} 
def get_nsamples (pool){
  return params.pools[pool].nsamples
}
def get_cr_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}
def get_cr_bam(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/possorted_genome_bam.bam", checkIfExists: true)
}


 
include { 
TEST_GZIP_INTEGRITY;
CELLRANGER;
CELLRANGER_VDJ;
SEURAT_PRE_FMX_QC
} from './pipeline_tasks.nf'

 
workflow  {
     
     // set up variables and channels
     library_dt = [] // channel with library and main data type (GEX or CITE)
     vdj_in = [] // channel of libraries with VDJ data
     no_vdj_in = [] // channel of libraries without VDJ data

     // read in the parameters
     list_pools = params.pools.keySet()
     for (pool in list_pools){
        libraries = get_libraries(pool)
        
        for (library in libraries){
          dts = get_data_types(pool, library)

          // check for CITE-seq data, if not the main data_type is "GEX"
          main_dt = ("CITE" in dts) ? "CITE" : "GEX"
          library_dt << [library, main_dt]
            
          // check for "TCR/BCR"
          if ("BCR" in dts){
            vdj_in << [library, "BCR"]
          } else {
            no_vdj_in << [library, "No BCR", "empty", "empty"]
          }
          if ("TCR" in dts){
            vdj_in << [library, "TCR"]
          } else {
            no_vdj_in << [library, "No TCR", "empty", "empty"]
          }

        }
     }
     
     
     // set up channels from input
     ch_all_lib_dt = Channel.fromList(library_dt)
      .mix(Channel.fromList(vdj_in)) // ['library', 'dt']
     ch_no_vdj = Channel.fromList(no_vdj_in)
 
    
     /* 
     --------------------------------------------------------
     Run cellranger
     --------------------------------------------------------
     */
     ch_vdj_libs = Channel.empty()
     ch_lib_h5 = Channel.empty()
     ch_library_dt = Channel.empty()
     if (!params.settings.skip_cellranger){ 

      TEST_GZIP_INTEGRITY(ch_all_lib_dt)
      ch_gzip_out = TEST_GZIP_INTEGRITY.out.branch{
        gex: it[1].contains("GEX")
        cite: it[1].contains("CITE")
        bcr: it[1].contains("BCR")
        tcr: it[1].contains("TCR")
      } // separate the output by data type
      
      // if it passes Gzip --> split the CITE/GEX to downstream channels
      ch_library_dt = ch_gzip_out.cite.mix(ch_gzip_out.gex) 
        .multiMap { it -> 
          cellranger_in: it
          prefilt_in: it
        } // [library, data_type]



      CELLRANGER(ch_library_dt.cellranger_in) // --> [[library, cr_bam, raw_h5], [cr_files]]
      ch_lib_h5 = CELLRANGER.out.bam_h5.map {it -> it[0,2]}
    
      // run cellranger for vdj if specified
      if (params.settings.add_tcr || params.settings.add_bcr ){
          CELLRANGER_VDJ(ch_gzip_out.bcr.mix(ch_gzip_out.tcr)) 
          
          // add "empty" TCR/BCR libs for missing
          ch_vdj_libs = CELLRANGER_VDJ.out.bam_h5
            .mix(ch_no_vdj)
            .branch { 
            tcr: it[1].contains("TCR")
            bcr: it[1].contains("BCR")
        }
      }
     } 
     // skip running cellranger and grab output from canonical location
     // TODO expand to work for VDJ as well
     else {
        lib_h5 = []
        list_pools = params.pools.keySet()
        for (pool in list_pools){
          for (library in get_libraries(pool)){
            h5 = get_cr_h5(library)
            lib_h5 << [library, h5]
          }
        }
        ch_lib_h5 = Channel.fromList(lib_h5)
        ch_library_dt = Channel.fromList(library_dt).multiMap { it -> prefilt_in: it}
     }

    
    /* 
    --------------------------------------------------------
    Filter barcode list for freemuxlet
    --------------------------------------------------------
    */

    ch_barcodes_list = Channel.empty()
    // filter before freemuxlet if specified
      ch_prefilt_in = ch_library_dt.prefilt_in
        .combine(ch_lib_h5, by:0)  // [library, data_type, raw_h5]
      SEURAT_PRE_FMX_QC(ch_prefilt_in) 

}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}







