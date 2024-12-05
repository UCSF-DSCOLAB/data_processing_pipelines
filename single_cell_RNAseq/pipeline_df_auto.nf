

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
SEURAT_ADD_BCR;
SEURAT_ADD_TCR;
SEURAT_LOAD_POST_QC;
SEURAT_POST_FILTER
} from './modules/pipeline_tasks.nf'


include {
get_c4_h5; get_library_ncells; get_libraries_data_type_tuples;
get_vdj_tuple; get_vdj_name; get_clonotypes; get_contigs; get_pre_fmx_qc_outputs;
get_pre_fmx_cutoffs; get_cutoffs; get_sample_map; get_sample_maps; get_c4_h5s
} from  './helpers/params_parse.nf'


include {
extractFileName
} from "./helpers/utils.nf"
 


workflow {


    ch_all_h5 = Channel.fromList(get_c4_h5s())
 
    ch_sample_map = Channel.fromList(get_sample_maps())

      /*
      --------------------------------------------------------
      Run doublet finder if specified
      --------------------------------------------------------
      */
     ch_initial_sobj = Channel.empty()
     if (params.settings.run_doubletfinder) {
        ch_doublet_input = Channel.from(get_library_ncells()).join(ch_all_h5).join(ch_sample_map) // [lib, ncells, raw_h5, fmx_cluster ]
        FIND_DOUBLETS(ch_doublet_input) // --> [lib, doublet_finder_sobj]
        ch_initial_sobj = FIND_DOUBLETS.out.sobj
     } else {
        ch_doublet_input = ch_all_h5.join(ch_sample_map) // [lib, ncells, raw_h5, fmx_cluster]
        LOAD_SOBJ(ch_doublet_input)
        ch_initial_sobj = LOAD_SOBJ.out.sobj
     }
     
      /* 
      --------------------------------------------------------
      Set up seurat object
      --------------------------------------------------------
      */
      // add TCR & BCR data
      
             ch_library_bcr_tcr = Channel.from(get_libraries_data_type_tuples()).transpose().filter { it[1] in ["BCR", "TCR"] }
                
              ch_vdj_libs = ch_library_bcr_tcr
              .map{
                it -> [it[0], it[1], get_clonotypes(it[0], it[1]), get_contigs(it[0], it[1])]
              }
              .branch { 
                        tcr: it[1].contains("TCR")
                        bcr: it[1].contains("BCR")
                    }  



     if (params.settings.add_tcr){
        SEURAT_ADD_TCR(ch_initial_sobj.join(ch_vdj_libs.tcr, by:0, remainder: true))
        ch_tcr_out = SEURAT_ADD_TCR.out.sobj
      } else {
        ch_tcr_out = ch_initial_sobj
      }

      if (params.settings.add_bcr){
        SEURAT_ADD_BCR(ch_tcr_out.join(ch_vdj_libs.bcr, by:0, remainder: true))
        ch_bcr_out = SEURAT_ADD_BCR.out.sobj
      } else {
        ch_bcr_out = ch_tcr_out
      } 

     /*
     --------------------------------------------------------
     Set up seurat object
     --------------------------------------------------------
     */
     ch_library_info = Channel.from(get_libraries_data_type_tuples()).transpose() // -> [[library_dir, data_type]]
     ch_seurat_input = ch_library_info.join(ch_bcr_out) // -> [library, data_type, ]
     .map{it -> [it[0], it[1], it[2], get_c4_h5(it[0]), get_pre_fmx_cutoffs(it[0])] } // pre_fmx
     SEURAT_LOAD_POST_QC(ch_seurat_input)

      // use cutoffs listed prior
      ch_seurat_post_qc_in = SEURAT_LOAD_POST_QC.out.cutoffs_file // -> [library, cutoffs]
          .combine(SEURAT_LOAD_POST_QC.out.qc_output, by:0) // -> [library, cutoffs, sobj
          .map{it -> [it[0], it[1], it[2], get_c4_h5(it[0])] } // -> [library, cutoffs, sobj, cr_h5]
      SEURAT_POST_FILTER(ch_seurat_post_qc_in)

}


