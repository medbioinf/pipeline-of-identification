nextflow.enable.dsl=2

// parameters for oktoberfest
params.oktoberfest_memory = "64 GB"
params.oktoberfest_intensity_model = "Prosit_2020_intensity_HCD"
params.oktoberfest_irt_mode = "Prosit_2019_irt"

/**
 * Runs oktoberfest rescoring for the given PSMs and mzML files.
 * 
 * @param psm_tsvs_and_mzmls: A tuple containing the PSM utils TSV files and the mzML files for the PSMs.
 * @param psm_tsvs: The PSM TSV files.
 * @param mzmls: The mzML files. 
 *
 * @return: The oktoberfest rescored PSMs in TSV format.
 */
workflow oktoberfest_rescore_workflow {
    take:
    psm_tsvs_and_mzmls
    psm_tsvs
    mzmls

    main:
    oktoberfest_pins = run_oktoberfest_feature_gen(psm_tsvs_and_mzmls, psm_tsvs, mzmls, params.fragment_tol_da)

    emit:
    oktoberfest_pins
}

/**
 * @param psm_tsvs_and_mzmls: A tuple containing the PSM utils TSV files and the mzML files for the PSMs.
 * @param psm_tsvs: The PSM TSV files.
 * @param mzmls: The mzML files. 
 * @param fragment_tol_da: The fragment tolerance for the rescoring.
 * 
 * @return The oktoberfest rescored PSMs in TSV format.
 */
process run_oktoberfest_feature_gen {
    cpus 1
    memory { params.oktoberfest_memory }

    container { params.python_image }

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val fragment_tol_da
    
    output:
    path "${psm_utils_tsvs.baseName}.features.tsv"
    
    script:
    """
    oktoberfest_feature_gen.py \
        -out-folder ./oktoberfest_out \
        -psms-file ${psm_utils_tsvs} \
        -spectra-file ${mzml_for_psms} \
        -intensity-model ${params.oktoberfest_intensity_model} \
        -irt-model ${params.oktoberfest_irt_modeö} \
        -mass-tolerance ${fragment_tol_da} \
        -mass-tolerance-unit da \
        -is-timstof ${params.is_timstof} \

    mv ./oktoberfest_out/results/none/rescore.tab "${psm_utils_tsvs.baseName}.features.tsv"

    // Clean up the output directory
    rm -r oktoberfest_out
    """
}
