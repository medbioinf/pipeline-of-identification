nextflow.enable.dsl=2


/**
 * Creates a default comet params file
 */
process create_default_comet_param_file {
    container 'medbioinf/ident-comparison-comet:latest'

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
    container 'medbioinf/ident-comparison-python:latest'

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
    maxForks 1
    container 'medbioinf/ident-comparison-comet:latest'

    input:
    val mzmls_and_comet_param_files
    path fasta
    path mzmls
    path comet_param_files

    output:
    path "*.txt"

    /** TODO: set the num_threads someway useful */
    """
    sed -i "s;^num_threads.*;num_threads = 16;" ${mzmls_and_comet_param_files[1]}

    sed -i "s;^output_sqtfile.*;output_sqtfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_txtfile.*;output_txtfile = 1;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_pepxmlfile.*;output_pepxmlfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_mzidentmlfile.*;output_mzidentmlfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^output_percolatorfile.*;output_percolatorfile = 0;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^print_expect_score.*;print_expect_score = 1;" ${mzmls_and_comet_param_files[1]}
    sed -i "s;^num_output_lines.*;num_output_lines = 5;" ${mzmls_and_comet_param_files[1]}
    
    comet -P${mzmls_and_comet_param_files[1]} -D${fasta} ${mzmls_and_comet_param_files[0]}
    """
}

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

        comet_tsvs = identification_with_comet(mzml_and_param_file, fasta, mzmls.collect(), comet_param_files.collect())
        comet_tsvs = comet_tsvs.flatten()
        
    emit:
        comet_tsvs
}