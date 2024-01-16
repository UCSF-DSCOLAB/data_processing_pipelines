def get_c4_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}

def get_c4_bam(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/possorted_genome_bam.bam", checkIfExists: true)
}

def get_cutoffs(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_cutoffs.csv", checkIfExists: true)
}

def get_sobj(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/automated_processing/${library}_raw.rds", checkIfExists: true)
}

def get_cr_h5(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cellranger/raw_feature_bc_matrix.h5", checkIfExists: true)
}

def get_c4_h5_bam(){
    return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_c4_bam(library.name), get_c4_h5(library.name)]
           }
    }
}

def get_pre_qc_outputs(){
	return params.pools.collectMany {
           pool -> pool.libraries.collect {
               library -> [library.name, get_cutoffs(library.name), get_sobj(library.name), get_cr_h5(library.name)]
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

def get_libraries_data_type(){
    return params.pools.collectMany {
                pool -> pool.libraries.collect {
                    library -> [library.name, library.data_types.join(",")]
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

