nextflow.enable.dsl=2

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (workflow.success){
       println "Deleting working directory $workDir"
       "rm -rf $workDir".execute()
    }
}

include { 
FIND_DOUBLETS;
LOAD_SOBJ;
SEURAT_QC;
SEURAT_POST_FILTER
} from './pipeline_tasks.nf'

def get_cutoffs(library){
    return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter/${library}_cutoffs.csv", checkIfExists: true)
}

def get_libraries (pool) {
  return params.pools[pool].libraries.keySet()
}
def get_data_types (pool, library){
  return params.pools[pool].libraries[library].data_types
} 

def get_cr_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}

def get_ncells (pool, library){
  return params.pools[pool].libraries[library].ncells_loaded
} 
def get_sample_map (library){
    if (params.settings.demux_method.equals("demuxlet")){   
        return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/demuxlet/${library}.clust1.samples.reduced.tsv", checkIfExists: true)
    } else {
        return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/freemuxlet/${library}.clust1.samples.reduced.tsv", checkIfExists: true)
    }
}


workflow {
   
    lib_ncells = []
    lib_h5 = []
    library_dt = []
    sample_map = []
    pool_libs = []
    list_pools = params.pools.keySet()
    for (pool in list_pools){
        libraries = get_libraries(pool)
        
        for (library in libraries){
          dts = get_data_types(pool, library)
          main_dt = ("CITE" in dts) ? "CITE" : "GEX"
          h5=get_cr_h5(library)
          pool_libs << [pool, library]
          sample_map << [library, get_sample_map(library)]
          lib_h5 << [library, h5]
          library_dt << [library, main_dt]
          ncells = get_ncells(pool, library)
          lib_ncells << [library, ncells]
        }
    }


    ch_lib_pools = Channel.fromList(pool_libs).map{ it-> it[1,0] }
    ch_lib_ncells = Channel.fromList(lib_ncells)
    ch_sample_map = Channel.fromList(sample_map)
    ch_library_dt = Channel.fromList(library_dt)


    // to set up:
    ch_lib_h5 = Channel.fromList(lib_h5)
    .multiMap { it -> 
        df_in: it 
        seurat_in: it
        post_seurat_in: it
    }
      
      /* 
      --------------------------------------------------------
      Run doublet finder if specified
      --------------------------------------------------------
      */
      ch_initial_sobj = Channel.empty()
      ch_df_in = ch_lib_h5.df_in.combine(ch_sample_map, by: 0) 
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
      // set up seurat QC
      ch_seurat_qc_in =  ch_library_dt //.seurat_in
      .combine(ch_initial_sobj, by:0)
      .combine(ch_lib_h5.seurat_in, by:0) // <-- [lib, data_type, doublet_finder_sobj, raw_h5]
      
      SEURAT_QC(ch_seurat_qc_in) 

        // use cutoffs listed prior
      lib_cutoffs = []
      for (pool in list_pools){
          for (library in get_libraries(pool)){
            cutoffs_file = get_cutoffs(library)
            lib_cutoffs << [library, cutoffs_file]
          }
        }

        ch_seurat_post_qc_in = Channel.fromList(lib_cutoffs)
          .combine(SEURAT_QC.out.qc_output, by:0)
          .combine(ch_lib_h5.post_seurat_in, by:0)
        SEURAT_POST_FILTER(ch_seurat_post_qc_in)

}


