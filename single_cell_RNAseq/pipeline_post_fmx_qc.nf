

nextflow.enable.dsl=2

include { 
SEURAT_PRE_FMX_FILTER;
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
SEURAT_QC;
FMX_ASSIGN_TO_GT;
SEURAT_POST_FILTER
} from './modules/pipeline_tasks.nf'

def get_cutoffs(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter/${library}_cutoffs.csv", checkIfExists: true)
}

def get_sobj(library){
  return file("${params.project_dir}/data/single_cell_GEX/processed/${library}/cell_filter/${library}_raw.rds", checkIfExists: true)
}

def get_libraries (pool) {
  return params.pools[pool].libraries.keySet()
}
def get_data_types (pool, library){
  return params.pools[pool].libraries[library].data_types
} 
def get_nsamples (pool){
  return params.pools[pool].nsamples
}
def get_vcf (pool){
  return params.pools[pool].vcf
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
def get_ncells (pool, library){
  return params.pools[pool].libraries[library].ncells_loaded
} 


workflow {
    library_dt = [] // channel with library and main data type (GEX or CITE)
    pool_libs = [] // list of all libraries
    lib_ncells = []
    pool_nsamples = [] // number of samples in a library (used for freemuxlet)
    pool_nlibs = [] // number of libraries per pool
    pool_vcf = []
    lib_h5 = []
    lib_bam = []
    post_qc_in = []
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
          h5=get_cr_h5(library)
          post_qc_in << [library, get_cutoffs(library), get_sobj(library)]
          main_dt = ("CITE" in dts) ? "CITE" : "GEX"
          library_dt << [library, main_dt]
          lib_h5 << [library, h5]
          lib_bam << [library, get_cr_bam(library)]

          ncells = get_ncells(pool, library)
          lib_ncells << [library, ncells]

        }
    }

    ch_pool_ns = Channel.fromList(pool_nsamples) // ['pool1', 5]
    ch_pool_libs = Channel.fromList(pool_libs)
    ch_lib_pools = Channel.fromList(pool_libs).map{ it-> it[1,0] }
    ch_pool_nlibs = Channel.fromList(pool_nlibs)
    ch_pool_vcf = Channel.fromList(pool_vcf)
    ch_library_dt = Channel.fromList(library_dt)
    ch_post_qc_in = Channel.fromList(post_qc_in)
    ch_lib_ncells = Channel.fromList(lib_ncells)
   

    // to set up:
    ch_lib_h5 = Channel.fromList(lib_h5)
    .multiMap { it -> 
        df_in: it 
        seurat_in: it
        post_seurat_in: it
    }
    ch_lib_bam = Channel.fromList(lib_bam)

    // set up if there's vdj
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
    // ch_cr_out.df_in [lib, raw_h5]
    // ch_cr_out.seurat_in


    SEURAT_PRE_FMX_FILTER(ch_post_qc_in)
    ch_barcodes_list = SEURAT_PRE_FMX_FILTER.out.bc_list


    // split the barcodes list into two channels for next steps
    ch_barcodes = ch_barcodes_list.multiMap{ it -> bam_in: dsc_in: it }


    /* 
    --------------------------------------------------------
    Setup bam, dsc files for freemuxlet
    --------------------------------------------------------
    */

    // Filter the bam file in prep for freemux
     FILTER_BAM(ch_lib_bam //ch_cr_out.filter_bam_in
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
          .groupTuple() // group by pool --> pool, libs, list_f
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
          pool_ns: it[0,1] // pool, nsamples
          pool_vcf: it[0,2] // pool, ref_vcf
        }
        ch_skip_merge_in_split = ch_merge_in.skip_merge.multiMap{ it ->
          lib_ns_files: it[0,1,3] // lib, nsamples
          lib_vcf: it[0,2] // lib, ref_vcf
        }
          
        MERGE_DSC(ch_merge_in_split.dsc_merge_in) // --> [pool1, pool1.tsv, merged_plp, merged_cl, merged_var, barcodes]       
        
        ch_merged = MERGE_DSC.out.merged_files.multiMap{ it -> 
          demux_in: it[0,2,3,4,5] // [pool1, merged_plp, merged_cl, merged_var, barcodes] 
          unmerge_in: it[0,1]  // [pool, pool_tsv]
        }


        if (params.settings.demux_method == "freemuxlet"){
          // run single library FMX on the single-lib ones
          FREEMUXLET_LIBRARY(ch_skip_merge_in_split.lib_ns_files)
          
          // run freemuxlet on merged
          ch_fmx_in = ch_merge_in_split.pool_ns
            .combine(ch_merged.demux_in, by: 0) // --> [pool1,  numsamples, ...]
          FREEMUXLET_POOL(ch_fmx_in) 

          if (params.settings.fmx_assign_to_gt){
            // TODO add checks that the reference file exists
              
            ch_ref_vcf = ch_skip_merge_in_split.lib_vcf.
              mix(ch_merge_in_split.pool_vcf) // [lib/pool, ref_vcf])
            ch_fmx_vcf = FREEMUXLET_LIBRARY.out.vcf.mix(FREEMUXLET_POOL.out.vcf) // [lib/pool, fmx_vcf]
            FMX_ASSIGN_TO_GT(ch_ref_vcf.combine(
              ch_fmx_vcf, by: 0) // [lib/pool, ref_vcf, fmx_vcf])
            )
          } 

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
        if (params.settings.demux_method == "demuxlet"){
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
         // ch_plp_files, ch_lib_pools, ch_pool_vcf
         if (params.settings.demux_method == "freemuxlet"){

          ch_fmx_in = ch_plp_files
              .combine(ch_lib_pools, by: 0)
              .map{it -> it[2,0,1]}
              .combine(ch_pool_ns, by:0) // add numsamples
              .map {it -> it[1,3,2]} //  --> [library, nsamples, files]
          FREEMUXLET_LIBRARY(ch_fmx_in)
          ch_sample_map = FREEMUXLET_LIBRARY.out.sample_map

          if (params.settings.fmx_assign_to_gt){
              // TODO add checks that the reference file exists
              FMX_ASSIGN_TO_GT(ch_pool_vcf.combine(
                FREEMUXLET_LIBRARY.out.vcf, by: 0) // [lib, ref_vcf, fmx_vcf])
              )
          } 


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
      ch_seurat_qc_in =  ch_library_dt //.seurat_in
      .combine(ch_bcr_out, by:0)
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


