nextflow.enable.dsl=2

msgfplus_image = 'quay.io/medbioinf/msgfplus:v2024.03.26'

// number of threads used by msamanda
params.msgfplus_threads = 16
params.msgfplus_mem_gb = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

/**
 * Executes the identification using MS-GF+
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow msgfplus_identification {
    take:
        msgfplus_params_file
        fasta
        mzmls
        precursor_tol_ppm

    main:
        msgfplus_results = identification_with_msgfplus(msgfplus_params_file, fasta, mzmls, precursor_tol_ppm)

        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msgfplus_results, 'mzid', 'msgfplus')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'msgfplus')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'msgfplus')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]
        
    emit:
        msgfplus_results
        psm_tsvs
        tda_results
        pin_files
        pout_files
        mokapot_results
        mokapot_features
}


process identification_with_msgfplus {
    cpus { params.msgfplus_threads }
    memory { params.msgfplus_mem_gb + " GB" }
    container { msgfplus_image }

    input:
    path msgfplus_params_file
    path fasta
    path mzmls
    val precursor_tol_ppm

    output:
    path "${mzmls.baseName}.mzid"

    """
    cp ${msgfplus_params_file} adjusted_MSGFPlus_Params.txt
    sed -i 's;PrecursorMassTolerance=.*;PrecursorMassTolerance=${precursor_tol_ppm};' adjusted_MSGFPlus_Params.txt

    java -Xmx${params.msgfplus_mem_gb}G -jar /opt/msgfplus/MSGFPlus.jar -conf adjusted_MSGFPlus_Params.txt -s ${mzmls} -d ${fasta} -thread ${params.msgfplus_threads} -o ${mzmls.baseName}.mzid
    """
}
