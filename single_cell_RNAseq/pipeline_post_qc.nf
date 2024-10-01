
nextflow.enable.dsl=2

workflow.onComplete {
//    println "Pipeline completed at: $workflow.complete"
//    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (workflow.success){
//      println "Deleting working directory $workDir"
       "rm -rf $workDir".execute()
    }
}


include { 
    SEURAT_POST_FILTER
} from './modules/pipeline_tasks.nf'

include {
get_pre_qc_outputs
} from  './helpers/params_parse.nf'


workflow {
    ch_pre_qc = Channel.fromList(get_pre_qc_outputs())
    SEURAT_POST_FILTER(ch_pre_qc) // [library, cutoffs, sobj, raw_h5]
}

