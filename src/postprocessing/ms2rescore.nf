nextflow.enable.dsl=2

// parameters for ms2rescore
params.ms2rescore_threads = 4
params.ms2rescore_mem = "64 GB"
params.ms2rescore_model = "HCD"
params.ms2rescore_chunk_size = 100000

workflow ms2rescore_workflow {
    take:
    psm_tsvs_and_mzmls
    psm_tsvs
    mzmls
    psm_id_pattern
    spectrum_id_pattern
    id_decoy_pattern
    searchengine

    main:
    ms2rescore_pre_pins = run_chunked_ms2rescore(psm_tsvs_and_mzmls, psm_tsvs, mzmls, spectrum_id_pattern, params.fragment_tol_da)
    ms2rescore_pins = correct_psm_utils_pins(ms2rescore_pre_pins)

    emit:
    ms2rescore_pins
}


process run_chunked_ms2rescore {
    cpus  { params.ms2rescore_threads }
    memory { params.ms2rescore_mem }

    container { params.python_image }
    containerOptions { "-v /mnt/data/projects/pipeline-of-identification/bin/ms2pip-model:/mnt/data/ms2pip-model" }

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val spectrum_id_pattern
    val fragment_tol_da
    
    output:
    path "*.ms2rescore.pin", emit: features_file
    
    script:
    """
    chunked_ms2rescore.py -psms_file ${psm_utils_tsvs} \
        -mzml_file ${mzml_for_psms} \
        -model ${params.ms2rescore_model} -model_dir "/mnt/data/ms2pip-model" \
        -ms2_tolerance ${fragment_tol_da} \
        -spectrum_id_pattern '${spectrum_id_pattern}' \
        -processes ${params.ms2rescore_threads} \
        -chunk_size ${params.ms2rescore_chunk_size} \
        -out_file "${psm_utils_tsvs}.ms2rescore.pin"
    """
}


process correct_psm_utils_pins {
    cpus  2
    memory '8 GB'

    container { params.python_image }

    input:
    path psm_utils_pins

    output:
    path "${psm_utils_pins.baseName}.corrected.pin"

    script:
    """
    # correct the PIN file by moving the scan number to third column and adding correct SpecId (increasing integer)
    awk '{FS="\t";OFS="\t"; if (NR>1) { \$3=\$1; \$1=NR-1; gsub(".*=", "", \$3)  } print}' ${psm_utils_pins} > ${psm_utils_pins.baseName}.corrected.pin
    """
}
