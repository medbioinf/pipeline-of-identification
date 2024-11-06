nextflow.enable.dsl=2

/**
 * This workflow uses the (enhanced) psm_utils TSV files as input and calculates the FDR q-values using the "default" TDA approach.
 */
workflow target_decoy_approach {
    take:
    psm_utils_tsvs
    searchengine

    main:
    filtered_small_dfs = calculate_qvalues_filter_psms_and_record_df(psm_utils_tsvs, searchengine)

    emit:
    filtered_small_dfs
}


/**
 * Caclulates the q-values and filters the PSMs (q-value < 0.01, only 1st spectrum rank identifications) based on the target-decoy approach.
 */
process calculate_qvalues_filter_psms_and_record_df {
    cpus 2
    memory '8GB'

    container { params.python_image }

    input:
    path psm_utils_tsv
    val searchengine

    output:
    path "${psm_utils_tsv.baseName}.tda.small_df.tsv"

    script:
    """
    tda_qvalues_and_filter_psm_tsvs.py -in_file ${psm_utils_tsv} -out_file ${psm_utils_tsv.baseName}.tda.small_df.tsv -searchengine ${searchengine}
    """
}
