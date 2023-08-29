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
         pools : ${params.pools.keySet()}
         """
         .stripIndent()


/* 
 * define getters for required tasks
 */

def get_libraries (pool) {
  return params.pools[pool].libraries.keySet()
}
def get_data_types (pool, library){
  return params.pools[pool].libraries[library].data_types
} 
def get_ncells (pool, library){
  return params.pools[pool].libraries[library].ncells_loaded
} 

def get_nsamples (pool){
  return params.pools[pool].nsamples
}
def get_vcf (pool){
  return params.pools[pool].vcf
}

def get_vdj_name(library, data_type){
  return library.replace("SCG", "SC" + data_type.substring(0, 1))
}
def get_vdj_tuple(library, data_type){
  return [library, data_type, get_vdj_name(library, data_type)]
}

def get_clonotypes(library, data_type){
  vdj_library=get_vdj_name(library, data_type)
  return file("${params.project_dir}/data/single_cell_${data_type}/processed/${vdj_library}/cellranger/clonotypes.csv")
}
def get_contigs(library, data_type){
  vdj_library=get_vdj_name(library, data_type)
  return file("${params.project_dir}/data/single_cell_${data_type}/processed/${vdj_library}/cellranger/filtered_contig_annotations.csv")
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

 
workflow  {
     
     // set up variables and channels
     pool_libs = [] // list of all libraries
     library_dt = [] // channel with library and main data type (GEX or CITE)
     vdj_in = [] // channel of libraries with VDJ data
     no_vdj_in = [] // channel of libraries without VDJ data
     pool_nsamples = [] // number of samples in a library (used for freemuxlet)
     pool_nlibs = [] // number of libraries per pool
     pool_vcf = []
     lib_ncells = []

     // read in the parameters
     list_pools = params.pools.keySet()
     for (pool in list_pools){
        libraries = get_libraries(pool)
        nsamples = get_nsamples(pool)
        vcf = get_vcf(pool)
        pool_nsamples << [pool, nsamples]
        pool_nlibs << [pool, libraries.size()]
        pool_vcf << [pool, vcf]
        
        for (library in libraries){
          dts = get_data_types(pool, library)
          pool_libs << [pool, library]
          // check for CITE-seq data, if not the main data_type is "GEX"
          main_dt = ("CITE" in dts) ? "CITE" : "GEX"
          library_dt << [library, main_dt]
          ncells = get_ncells(pool, library)
          lib_ncells << [library, ncells]
            
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
     ch_pool_ns = Channel.fromList(pool_nsamples) // ['pool1', 5]
     ch_pool_libs = Channel.fromList(pool_libs)
     ch_lib_pools = Channel.fromList(pool_libs).map{it-> it[1,0]}
     ch_pool_nlibs = Channel.fromList(pool_nlibs)
     ch_no_vdj = Channel.fromList(no_vdj_in)
     ch_pool_vcf = Channel.fromList(pool_vcf)
     ch_lib_ncells = Channel.fromList(lib_ncells)
 
    
     /* 
     --------------------------------------------------------
     Run cellranger
     --------------------------------------------------------
     */
     ch_cr_bam_h5 = Channel.empty()
     ch_vdj_libs = Channel.empty()

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
          seurat_in: it 
        } // [library, data_type]



      CELLRANGER(ch_library_dt.cellranger_in) // --> [[library, cr_bam, raw_h5], [cr_files]]
      ch_cr_bam_h5 = CELLRANGER.out.bam_h5
    
      // run cellranger for vdj if specified
      if (params.settings.add_tcr || params.settings.add_bcr ){
          ch_vdj_in = ch_gzip_out.bcr.mix(ch_gzip_out.tcr).map{
            it -> get_vdj_tuple(it[0], it[1])
          }
          CELLRANGER_VDJ(ch_vdj_in) 
          
          // add "empty" TCR/BCR libs for missing
          ch_vdj_libs = CELLRANGER_VDJ.out.vdj_csvs
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

        cr_bam_h5 = []
        cr_vdj = []
        list_pools = params.pools.keySet()
        for (pool in list_pools){
          for (library in get_libraries(pool)){
            bam = get_cr_bam(library)
            h5 = get_cr_h5(library)
            cr_bam_h5 << [library, bam, h5]
            

          }
        }

        ch_cr_bam_h5 = Channel.fromList(cr_bam_h5)
        ch_library_dt = Channel.fromList(library_dt).multiMap { it -> seurat_in: it}


        // vdj_in [library, dt]

        if (params.settings.add_tcr || params.settings.add_bcr ){
            ch_vdj_in = Channel.fromList(vdj_in).map{
              it -> [it[0], it[1], get_clonotypes(it[0], it[1]), get_contigs(it[0], it[1])]
            }

            ch_vdj_libs = ch_vdj_in.mix(ch_no_vdj)
                    .branch { 
                        tcr: it[1].contains("TCR")
                        bcr: it[1].contains("BCR")
                    }  
        }
     }

     ch_cr_out = ch_cr_bam_h5 
      .multiMap{ it ->
        barcodes_in: it[0,2] //  [library, raw_h5]
        filter_bam_in: it[0,1]  // [library, cr_bam]
        df_in: it[0,2] //  [library, raw_h5]
        seurat_in: it[0,2] //  [library, raw_h5]
      }



    
    /* 
    --------------------------------------------------------
    Filter barcode list for freemuxlet
    --------------------------------------------------------
    */

    ch_barcodes_list = Channel.empty()

    FILTER_BARCODES(ch_cr_out.barcodes_in) 
    ch_barcodes_list = FILTER_BARCODES.out.bc_list // [library, barcodes]
    
    
     // split the barcodes list into two channels for next steps
     ch_barcodes = ch_barcodes_list.multiMap{ it -> bam_in: dsc_in: it }


    /* 
    --------------------------------------------------------
    Setup bam, dsc files for freemuxlet
    --------------------------------------------------------
    */

    // Filter the bam file in prep for freemux
     FILTER_BAM(ch_cr_out.filter_bam_in
      .join(ch_barcodes.bam_in) 
      ) // [library, cr_bam, barcodes] --> [library, filtered_bam]
     
     // run dsc_pileup
     DSC_PILEUP(ch_barcodes.dsc_in
      .join(FILTER_BAM.out.bam_file)  // [library, barcodes, filtered_bam]
      )
     ch_plp_files = DSC_PILEUP.out.plp_files
     /* 
     --------------------------------------------------------
     Run freemuxlet
     --------------------------------------------------------
     */

     ch_sample_map = Channel.empty() 

     // ---- merge freemuxlet by pool if specified ---- //
     if (params.settings.merge_for_demux){
        
        // prepare input channel: first group by pool
        ch_dsc_by_pool = ch_plp_files
          .combine(ch_lib_pools, by:0) // add pool
          .map { it ->  it[2,0,1] } // pool, lib, list_f
          .groupTuple() // group by pool 
          .combine(ch_pool_nlibs, by:0) // add nlibs
          .combine(ch_pool_ns, by:0) // add num samples for freemux
          .combine(ch_pool_vcf, by:0) // add vcf for demux --> [pool, libs, list_f, nlibs, nsamples, vcf]
          .map {it -> it[0,3,1,4,5,2]} // --> [pool, nlibs, libs, nsamples, vcf, list_f]
                  
        // second separate out the ones to merge (multi-lib pool) & not (single lib pool)
        ch_merge_in = ch_dsc_by_pool
          .branch{
            merge: it[1] > 1 // multi-lib
              return tuple(it[0], it[3], it[4], it[5].flatten().collect()) // pool, nsamples, vcf, files
            skip_merge: it[1] == 1 // single lib
              return tuple(it[2][0], it[3], it[4], it[5][0]) // libs, nsamples, vcf, files
         }
        

        // run merge on the multi-lib ones 
        // first separate for dsc vs fmx steps
        ch_merge_in_split = ch_merge_in.merge.multiMap{ it ->
          dsc_merge_in: it[0,3] // pool, files
          pool_ns_vcf: it[0,1,2] // pool, nsamples, vcf
        }
          
        MERGE_DSC(ch_merge_in_split.dsc_merge_in) // --> [pool1, pool1.tsv, merged_plp, merged_cl, merged_var, barcodes]       
        
        ch_merged = MERGE_DSC.out.merged_files.multiMap{ it -> 
          demux_in: it[0,2,3,4,5] // [pool1, merged_plp, merged_cl, merged_var, barcodes] 
          unmerge_in: it[0,1]  // [pool, pool_tsv]
        }
        
        if (params.settings.demux_method == "freemuxlet"){
          // run single library FMX on the single-lib ones
          FREEMUXLET_LIBRARY(ch_merge_in.skip_merge.map{
            it -> it[0,1,3] // remove VCF
            })
          
          // run freemuxlet on merged
          ch_fmx_in = ch_merge_in_split.pool_ns_vcf
            .map{ it -> it[0,1]} // remove vcf
            .combine(ch_merged.demux_in, by: 0) // --> [pool1,  numsamples, ...]
          FREEMUXLET_POOL(ch_fmx_in) 

          // separate demulitplexing output files
          ch_unmerge_in = ch_merged.unmerge_in.combine( FREEMUXLET_POOL.out.merged_files, by: 0)
          UNMERGE_FMX(ch_unmerge_in) // --> [pool, list_of_f]
          
          ch_separate = UNMERGE_FMX.out.samples_file
            .flatten()
            .map { it -> 
            def lib = it.name.toString().tokenize('.').get(0) // extracts the library
            return tuple(lib, it)
            }
            .groupTuple() // group by library

          SEPARATE_FMX(ch_separate) // --> lib, fmx_clus
          ch_sample_map = SEPARATE_FMX.out.sample_map
            .mix(FREEMUXLET_LIBRARY.out.sample_map)
        }
        if (params.settings.demux_method == "demuxlet") {
          DEMUXLET_LIBRARY(ch_merge_in.skip_merge.map{
            it -> it[0,2,3] // remove nsamples
            })

          // run freemuxlet on merged
          ch_dmx_in = ch_merge_in_split.pool_ns_vcf
            .map{ it -> it[0,2]} // remove nsamples
            .combine(ch_merged.demux_in, by: 0) // --> [pool1,  vcf, ...]
          DEMUXLET_POOL(ch_dmx_in) 

          ch_separate = DEMUXLET_POOL.out.merged_best
            .combine(ch_pool_libs, by: 0)
            .map{ it[2,1] } // get rid of pool

          SEPARATE_DMX(ch_separate)
          ch_sample_map = SEPARATE_DMX.out.sample_map
            .mix(DEMUXLET_LIBRARY.out.sample_map)
        }
      } 

       // ---- if not merging, run on each library individually --- //
       else { 
         if (params.settings.demux_method == "freemuxlet"){
          ch_fmx_in = ch_plp_files
              .combine(ch_lib_pools, by: 0)
              .map{it -> it[2,0,1]}
              .combine(ch_pool_ns, by:0) // add numsamples
              .map {it -> it[1,3,2]} //  --> [library, nsamples, files]
          FREEMUXLET_LIBRARY(ch_fmx_in)
          ch_sample_map = FREEMUXLET_LIBRARY.out.sample_map
         }
         if (params.settings.demux_method == "demuxlet"){
          ch_dmx_in = ch_plp_files
              .combine(ch_lib_pools, by: 0)
              .map{it -> it[2,0,1]}
              .combine(ch_pool_vcf, by:0) // add vcf
              .map {it -> it[1,3,2]} //  --> [library, vcf, files]
            
          DEMUXLET_LIBRARY(ch_dmx_in)
          ch_sample_map = DEMUXLET_LIBRARY.out.sample_map
         }

      }
      
      /* 
      --------------------------------------------------------
      Run doublet finder if specified
      --------------------------------------------------------
      */
      ch_initial_sobj = Channel.empty()
      ch_df_in = ch_cr_out.df_in.combine(ch_sample_map, by: 0) // --> [lib, ncells, raw_h5, fmx_clusters]
      if (params.settings.run_doubletfinder) {
        ch_df_in2 = ch_lib_ncells.combine(ch_df_in, by: 0) 
        FIND_DOUBLETS(ch_df_in2) // --> [lib, doublet_finder_sobj]
        ch_initial_sobj = FIND_DOUBLETS.out.sobj
      } 
      else { 
        LOAD_SOBJ(ch_df_in) 
        ch_initial_sobj = LOAD_SOBJ.out.sobj
      }

            /* 
      --------------------------------------------------------
      Set up seurat object
      --------------------------------------------------------
      */
      // add TCR & BCR data
      ch_tcr_out = Channel.empty()
      ch_bcr_out = Channel.empty()
      if (params.settings.add_tcr){
        SEURAT_ADD_TCR(ch_initial_sobj.combine(ch_vdj_libs.tcr, by:0))
        ch_tcr_out = SEURAT_ADD_TCR.out.sobj
      } else {
        ch_tcr_out = ch_initial_sobj
      }

      if (params.settings.add_bcr){
        SEURAT_ADD_BCR(ch_tcr_out.combine(ch_vdj_libs.bcr, by:0))
        ch_bcr_out = SEURAT_ADD_BCR.out.sobj
      } else {
        ch_bcr_out = ch_tcr_out
      }
      
      // set up seurat QC
      ch_seurat_qc_in = ch_library_dt.seurat_in
      .combine(ch_bcr_out, by:0)
      .combine(ch_cr_out.seurat_in, by:0) // <-- [lib, data_type, doublet_finder_sobj, raw_h5]
      
      SEURAT_QC(ch_seurat_qc_in) 
      
  
}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}







