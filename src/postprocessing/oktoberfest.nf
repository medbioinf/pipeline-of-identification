nextflow.enable.dsl=2

// parameters for oktoberfest
params.oktoberfest_memory = "64 GB"
params.oktoberfest_intensity_model = "Prosit_2020_intensity_HCD"
params.oktoberfest_irt_model = "Prosit_2019_irt"

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
    oktoberfest_features = run_oktoberfest_feature_gen(psm_tsvs_and_mzmls, psm_tsvs, mzmls, params.fragment_tol_da)
    oktoberfest_pins = oktoberfest_features_to_pin(oktoberfest_features)


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
        -irt-model ${params.oktoberfest_irt_model} \
        -mass-tolerance ${fragment_tol_da} \
        -mass-tolerance-unit da \
        -is-timstof ${params.is_timstof} \

    mv ./oktoberfest_out/results/none/rescore.tab "${psm_utils_tsvs.baseName}.features.tsv"

    // Clean up the output directory
    rm -r oktoberfest_out
    """
}

/**
 * @param okt_features_tsv: Oktoberfest feature file.
 * 
 * @return Oktoberfest feature file in PIN format ready to use with percolator.
 */
process oktoberfest_features_to_pin {
    cpus 1
    memory { params.oktoberfest_memory }

    container { params.python_image }

    input:
    path okt_features_tsv

    output:
    path "${okt_features_tsv.baseName}.oktoberfest.pin"

    script:
    """
    oktoberfest_feature_to_pin.py \
        -features-file ${okt_features_tsv} \
        -out-folder ./{okt_features_tsv.baseName}.oktoberfest.pin
    """
}
