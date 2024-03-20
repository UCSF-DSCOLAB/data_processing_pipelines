nextflow.enable.dsl=2

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (workflow.success){
       println "Deleting working directory $workDir"
       "rm -rf $workDir".execute()
    }
}


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


 
include { 
TEST_GZIP_INTEGRITY;
CELLRANGER;
CELLRANGER_VDJ;
SEURAT_PRE_FMX_QC
} from './modules/pipeline_tasks.nf'

include {
get_c4_h5_bam; get_pool_library_meta; get_libraries_data_type_tuples;
} from  './helpers/params_parse.nf'


include {
extractFileName
} from "./helpers/utils.nf"
 
workflow  {
     
      ch_gex_cite_bam_h5 = Channel.empty()

      if (params.settings.skip_cellranger){
            ch_gex_cite_bam_h5 =  Channel.from(get_c4_h5_bam()) // [[library, cell_ranger_bam, raw_h5]
            
      } else {
            ch_library_info = Channel.from(get_libraries_data_type_tuples()).transpose()
            TEST_GZIP_INTEGRITY(ch_library_info) // -> [[library_dir, data_type]]
            
            ch_gzip_out = TEST_GZIP_INTEGRITY.out
            .branch{
              gex_cite: it[1] in ["GEX", "CITE"]
              bcr_tcr: it[1] in ["BCR", "TCR"]
            }

            // Run cellranger for GEX and CITE data types
            CELLRANGER(ch_gzip_out.gex_cite)
            ch_gex_cite_bam_h5 = CELLRANGER.out.bam_h5 // --> [[library, cell_ranger_bam, raw_h5]]

            // Run cellranger for BCR and TCR data types
            if (params.settings.add_tcr || params.settings.add_bcr ){
                
                ch_vdj_in = ch_gzip_out.bcr_tcr.map{
                  it -> get_vdj_tuple(it[0], it[1])
                }
                CELLRANGER_VDJ(ch_vdj_in) 
 
            }
      }

        ch_all_h5 = ch_gex_cite_bam_h5.map { it -> [it[0], it[2]] } // [[library, raw_h5 ]]

        ch_main_dt = Channel.from(get_libraries_data_type_tuples()).transpose()
            .filter(it -> it[1] in ["GEX", "CITE"])

        ch_prefilt_in = ch_main_dt.join(ch_all_h5, by:0) // [library, data_type, h5]

        SEURAT_PRE_FMX_QC(ch_prefilt_in) 

}








