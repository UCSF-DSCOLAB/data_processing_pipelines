def get_c4_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}

def get_c4_cr_filt_bc(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/filtered_feature_bc_matrix/barcodes.tsv.gz", checkIfExists: true)
}

def get_c4_bam(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/possorted_genome_bam.bam", checkIfExists: true)
}

def get_cutoffs(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_cutoffs.csv", checkIfExists: true)
}

def get_pre_fmx_cutoffs(library){
  path = "${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter/${library}_cutoffs.csv"
  if (file(path).exists()){
    return file(path, checkIfExists: true)
  } else {
    return(get_cutoffs(library))
  }
}

def get_sobj(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_raw.rds", checkIfExists: true)
}

def get_pre_fmx_sobj(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter/${library}_raw.rds", checkIfExists: true)
}


def get_c4_h5_bam(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_c4_bam(library.name), get_c4_h5(library.name)]
           }
    }
}

def get_c4_h5s(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_c4_h5(library.name)]
    } 
  }
}

def get_c4_h5_bam_bc(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_c4_bam(library.name), get_c4_h5(library.name), get_c4_cr_filt_bc(library.name)]
           }
    }
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
  return file("${params.project_dir}/data/single_cell_${data_type}/processed/${vdj_library}/cellranger/all_contig_annotations.csv")
}
def get_sample_map(library){
    return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/${params.settings.demux_method}/${library}.clust1.samples.reduced.tsv")
}

def get_sample_maps(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_sample_map(library.name)]
           }
    }
}



def get_pre_qc_outputs(){
	return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_cutoffs(library.name), get_sobj(library.name), get_c4_h5(library.name)]
           }
    }
}

def get_pre_fmx_qc_outputs(){
	return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_pre_fmx_cutoffs(library.name), get_pre_fmx_sobj(library.name), get_c4_h5(library.name)]
           }
    }
}


def get_pool_library_meta(){
    return params.pools.collectMany {
        pool ->
        return [
            [
                name: pool.name,
                vcf: pool.vcf,
                nsamples: pool.nsamples,
                num_of_libraries: pool.libraries.size(),
                lib_directories: pool.libraries*.name
            ]
        ]
    }
}



def get_libraries_data_type_tuples(){
    return params.pools.collectMany {
                pool -> pool.libraries.collect {
                    library -> [library.name, library.data_types]
                }
        }
           
}


def get_library_ncells(){
    return params.pools.collectMany {
                pool -> pool.libraries.collect {
                    library -> [library.name, library.ncells_loaded]
                }
           }
}

def get_multi_library_by_pool() {
    return get_pool_library_meta().findAll {it.num_of_libraries > 1}.collectMany { pool ->
            pool.lib_directories.collect { name ->
                [name, pool.name]
            }
        }
}

def get_multi_pool_by_library() {
    return get_pool_library_meta().findAll {it.num_of_libraries > 1}.collectMany { pool ->
            pool.lib_directories.collect { name ->
                [pool.name, name]
            }
        }
}



def get_single_library_by_pool() {
    return get_pool_library_meta().findAll {it.num_of_libraries == 1}.collectMany { pool ->
            pool.lib_directories.collect { name ->
                [name, pool.name]
            }
        }
}

def get_library_by_pool() {
    return get_pool_library_meta().collectMany { pool ->
            pool.lib_directories.collect { name ->
                [name, pool.name]
            }
        }
}


def get_library_by_sample_count() {
    return get_pool_library_meta().collectMany { pool ->
            pool.lib_directories.collect { name ->
                [name, pool.nsamples]
            }
        }
}

def get_pool_by_sample_count() {
    return get_pool_library_meta().collectMany { pool ->
            [
                [pool.name, pool.nsamples]
            ]
        }
}

def get_pool_vcf() {
    return get_pool_library_meta().collectMany { pool ->
            [
                [pool.name, pool.vcf]
            ]
        }

}

