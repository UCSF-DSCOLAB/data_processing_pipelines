
nextflow.enable.dsl=2

include { 
    SEURAT_POST_FILTER
} from './pipeline_tasks.nf'

def get_cutoffs(library){
 if (params.settings.demux_method.equals("demuxlet")){
        return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing_dmx/${library}_cutoffs.csv", checkIfExists: true)
    } else {
        return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_cutoffs.csv", checkIfExists: true)
    }
}

def get_sobj(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_raw.rds", checkIfExists: true)
}

def get_cr_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}

def get_libraries (pool) {
  return params.pools[pool].libraries.keySet()
}

workflow {

    post_qc_in = []
    list_pools = params.pools.keySet()
    for (pool in list_pools){
        libraries = get_libraries(pool)   
        for (library in libraries){
          post_qc_in << [library, get_cutoffs(library), get_sobj(library), get_cr_h5(library)]
        }
    }
    ch_post_qc_in = Channel.fromList(post_qc_in)
    SEURAT_POST_FILTER(ch_post_qc_in) // [library, cutoffs, sobj, raw_h5]
}

