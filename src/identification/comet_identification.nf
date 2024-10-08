nextflow.enable.dsl=2

comet_image = 'quay.io/medbioinf/comet-ms:v2024.01.0'
python_image = 'medbioinf/ident-comparison-python'

// number of threads used by comet
params.comet_threads = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

/**
 * Exports the identification using Comet configured by a SDRF files
 */
workflow comet_identification {
    take:
        sdrf
        fasta
        mzmls

    main:
        default_comet_params_file = create_default_comet_param_file()

        comet_param_files = create_comet_param_files_from_sdrf(sdrf, default_comet_params_file)
        comet_param_files = comet_param_files.flatten()
        
        // create a map containing the mzML file and the corresponding comet_param_file
        mzml_and_param_file = comet_param_files.map { it -> [ it.name.take(it.name.lastIndexOf('.comet.params')), it.name ] }

        comet_mzids = identification_with_comet(mzml_and_param_file, fasta, mzmls.collect(), comet_param_files.collect())
        comet_mzids = comet_mzids.flatten()
        
        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(comet_mzids, 'mzid', 'comet')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'comet')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'comet')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]

    emit:
        comet_mzids
        psm_tsvs
        tda_results
        pin_files
        pout_files
        mokapot_results
        mokapot_features
}

/**
 * Creates a default comet params file
 */
process create_default_comet_param_file {
    container { comet_image }

    output:
    path "comet.params.new"

    """
    comet -p
    """
}

/**
 * Creates a comet params file from an SDRF file
 * @param sdrf The SDRF file
 * @param default_comet_params_file The default comet params file
 * @return The comet params file
 */
process create_comet_param_files_from_sdrf {
    container { python_image }

    input:
    path sdrf
    path default_comet_params_file

    output:
    path "*.params"

    """
    python -m sdrf_convert $sdrf comet --config-folder ./ $default_comet_params_file
    """
}


process identification_with_comet {
    cpus { params.comet_threads }
    container { comet_image }

    input:
    val mzmls_and_comet_param_files
    path fasta
    path mzmls
    path comet_param_files

    output:
    path "*.mzid"

    """
    sed -i "s;^num_threads.*;num_threads = ${params.comet_threads};" ${mzmls_and_comet_param_files[1]}

    sed -i "s;^output_sqtfile.*;output_sqtfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_txtfile.*;output_txtfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_pepxmlfile.*;output_pepxmlfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_mzidentmlfile.*;output_mzidentmlfile = 1;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_percolatorfile.*;output_percolatorfile = 0;" ${mzmls_and_comet_param_files[1]}
    
    sed -i "s;^num_output_lines.*;num_output_lines = 5;" ${mzmls_and_comet_param_files[1]}
    
    comet -P${mzmls_and_comet_param_files[1]} -D${fasta} ${mzmls_and_comet_param_files[0]}
    """
}
