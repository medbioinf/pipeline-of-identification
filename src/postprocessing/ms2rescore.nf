nextflow.enable.dsl=2

ms2rescore_image = 'medbioinf/ident-comparison-ms2rescore'

// number of threads used by ms2rescore
params.ms2rescore_threads = 4
params.ms2rescore_mem = "64 GB"


workflow ms2rescore_workflow {
    take:
        psm_tsvs_and_mzmls
        psm_tsvs
        mzmls
        searchengine

    main:
        ms2rescore_results = run_ms2rescore(psm_tsvs_and_mzmls, psm_tsvs, mzmls, searchengine)

    emit:
        ms2rescore_results.mokapot_results
        ms2rescore_results.features_file
}


process run_ms2rescore {
    cpus  { params.ms2rescore_threads }
    memory { params.ms2rescore_mem }

    container { ms2rescore_image }
    containerOptions { "-v /mnt/data/projects/pipeline-of-identification/bin/ms2pip-model:/mnt/data/ms2pip-model" }

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val searchengine

    output:
    path "*.ms2rescore.mokapot.psms.txt", emit: mokapot_results
    path "*.ms2rescore.psms.tsv", emit: features_file

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
model = "HCD"
ms2_tolerance = 0.02
processes = ${params.ms2rescore_threads}
model_dir = "/mnt/data/ms2pip-model"

[ms2rescore.rescoring_engine.mokapot]
write_weights = true
write_txt = true
write_flashlfq = false' >> ms2rescore-config.toml

    ms2rescore -c ms2rescore-config.toml
    """
}
