nextflow.enable.dsl=2

/**
 * Creates a default comet params file
 */
process create_default_comet_params_file {
    container 'medbioinforub/ident-comparison-comet:latest'

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
process create_comet_params_file_from_sdrf {
    container 'medbioinforub/ident-comparison-python:latest'

    input:
    path sdrf
    path default_comet_params_file

    output:
    path "*.params"

    """
    python -m sdrf_convert $sdrf comet <MAX_VAFRIABLE_MODS> $default_comet_params_file .
    """
}

process identification_with_comet {
    maxForks 1
    container 'medbioinforub/ident-comparison-comet:latest'

    input:
    path comet_param_file
    path fasta
    path mzml

    output:
    path "*.mzIdentML"

    """
    comet -P${comet_param_file} -D${fasta} ${mzml}
    """
}

/**
 * Exports the identification using Comet configured by a SDRF files
 */
workflow convert_raw_files_to_mzml {
    take:
        mzml

    main:
        default_comet_params_file = create_default_comet_params_file()
        comet_param_files = create_comet_params_file_from_sdrf(sdrf)
        mzidentmls = identification_with_comet(thermo_raw_files)
    emit:
        mzidentmls
}