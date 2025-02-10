nextflow.enable.dsl=2

params.ms2rescore_image = 'medbioinf/ident-comparison-ms2rescore'

// parameters for ms2rescore
params.ms2rescore_threads = 4
params.ms2rescore_mem = "64 GB"
params.ms2rescore_model = "HCD"

// include filtering from convert_and_enhance_psm_tsv.nf
include {filter_pin_keep_only_best} from './convert_and_enhance_psm_tsv.nf'

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
    ms2rescore_pre_pins = run_ms2rescore(psm_tsvs_and_mzmls, psm_tsvs, mzmls, psm_id_pattern, spectrum_id_pattern, id_decoy_pattern, params.fragment_tol_da)
    ms2rescore_pins = correct_psm_utils_pins(ms2rescore_pre_pins)
    ms2rescore_onlybest_base_pins = filter_pin_keep_only_best(ms2rescore_pins, searchengine)

    emit:
    ms2rescore_pins
    ms2rescore_onlybest_base_pins
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
    val psm_id_pattern
    val spectrum_id_pattern
    val id_decoy_pattern
    val fragment_tol_da

    output:
    path "*.ms2rescore.pin", emit: features_file

    script:
    """
    # write the toml file...
    echo "[ms2rescore]
psm_file = '${psm_utils_tsvs}'
psm_file_type = 'tsv'
spectrum_path = '${mzml_for_psms}'

psm_id_pattern = '${psm_id_pattern}'             # scan/spectrum id in PSM file
spectrum_id_pattern = '${spectrum_id_pattern}'   # scan/spectrum id in mzML, Single quotes for literal regex string" > ms2rescore-config.toml

    # don't use DECOY_ prefix for MaxQuant, but all other searchengines
    if [ "" != "${id_decoy_pattern}" ]; then
        echo "id_decoy_pattern = '${id_decoy_pattern}'" >> ms2rescore-config.toml
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
    awk '{FS="\t";OFS="\t"; if (NR>1) { \$3=\$1; \$1=NR-1; gsub(".*=", "", \$3)  } print}' ${psm_utils_pins} > ${psm_utils_pins.baseName}.corrected.pin
    """
}
