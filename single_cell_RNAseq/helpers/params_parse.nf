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

def get_pool_library_meta(){
    return params.pools.collectMany {
        pool ->
        return [
            [
                name: pool.name,
                nsamples: pool.nsamples,
                num_of_libraries: pool.libraries.size(),
                lib_directories: pool.libraries*.dir
            ]
        ]
    }
}

def get_libraries_data_type(){
    return params.pools.collectMany {
                pool -> pool.libraries.collect {
                    library -> [library.dir, library.data_types.join(",")]
                }
           }
}

def get_pools_with_multi_library() {
    return get_pool_library_meta().findAll {it.num_of_libraries > 1}.collectMany { pool ->
            pool.lib_directories.collect { dir ->
                pool.name
            }
        }.unique()
}

def get_multi_library_by_pool() {
    return get_pool_library_meta().findAll {it.num_of_libraries > 1}.collectMany { pool ->
            pool.lib_directories.collect { dir ->
                [dir, pool.name]
            }
        }
}

def get_single_library_by_pool() {
    return get_pool_library_meta().findAll {it.num_of_libraries == 1}.collectMany { pool ->
            pool.lib_directories.collect { dir ->
                [dir, pool.name]
            }
        }
}


def get_library_by_sample_count() {
    return get_pool_library_meta().collectMany { pool ->
            pool.lib_directories.collect { dir ->
                [dir, pool.nsamples]
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


