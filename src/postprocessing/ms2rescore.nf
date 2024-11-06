nextflow.enable.dsl=2

params.ms2rescore_image = 'medbioinf/ident-comparison-ms2rescore'

// number of threads used by ms2rescore
params.ms2rescore_threads = 4
params.ms2rescore_mem = "64 GB"

params.ms2rescore_model = "HCD"

workflow ms2rescore_workflow {
    take:
    psm_tsvs_and_mzmls
    psm_tsvs
    mzmls
    searchengine

    main:
    ms2rescore_pre_pins = run_ms2rescore(psm_tsvs_and_mzmls, psm_tsvs, mzmls, searchengine, params.fragment_tol_da)
    ms2rescore_pins = correct_psm_utils_pins(ms2rescore_pre_pins)

    emit:
    ms2rescore_pins
}


process run_ms2rescore {
    cpus  { params.ms2rescore_threads }
    memory { params.ms2rescore_mem }

    container { params.ms2rescore_image }
    containerOptions { "-v /mnt/data/projects/pipeline-of-identification/bin/ms2pip-model:/mnt/data/ms2pip-model" }

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val searchengine
    val fragment_tol_da

    output:
    path "*.ms2rescore.pin", emit: features_file

    script:
    """
    # write the toml file...
    echo '[ms2rescore]
psm_file = "${psm_utils_tsvs}"
psm_file_type = "tsv"
spectrum_path = "${mzml_for_psms}"

psm_id_pattern = "(.*)"                 # scan/spectrum id in PSM file' > ms2rescore-config.toml


    if [ "${searchengine}" == "comet" ] || [ "${searchengine}" == "maxquant" ]; then
        echo "spectrum_id_pattern = '.*scan=(\\d+)\$'   # scan/spectrum id in mzML, Single quotes for literal regex string" >> ms2rescore-config.toml
    else
        echo "spectrum_id_pattern = '(.*)'   # scan/spectrum id in mzML, Single quotes for literal regex string" >> ms2rescore-config.toml
    fi

    # don't use DECOY_ prefix for MaxQuant, but all other searchengines
    if [ "${searchengine}" != "maxquant" ]; then
        echo "id_decoy_pattern = '^DECOY_'" >> ms2rescore-config.toml
    fi

    echo 'lower_score_is_better = false

max_psm_rank_input = 5
max_psm_rank_output = 1

processes = ${params.ms2rescore_threads}

[ms2rescore.feature_generators.ms2pip]
model = "${params.ms2rescore_model}"
ms2_tolerance = ${fragment_tol_da}
processes = ${params.ms2rescore_threads}
model_dir = "/mnt/data/ms2pip-model"

[ms2rescore.rescoring_engine]
# empty, so no rescoring with percolator or mokapot is performed
' >> ms2rescore-config.toml

    ms2rescore -c ms2rescore-config.toml
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
    awk '{FS="\t";OFS="\t"; if (NR>1) { \$3=\$1; \$1=NR-1; gsub(".*scan=", "", \$3)  } print}' ${psm_utils_pins} > ${psm_utils_pins.baseName}.corrected.pin
    """
}
