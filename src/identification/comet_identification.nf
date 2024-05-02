nextflow.enable.dsl=2


/**
 * Creates a default comet params file
 */
process create_default_comet_params_file {
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
process create_comet_params_file_from_sdrf {
    container 'medbioinf/ident-comparison-python:latest'

    input:
    path sdrf
    path default_comet_params_file
    val max_variable_mods

    output:
    path "*.params"

    """
    python -m sdrf_convert $sdrf comet $max_variable_mods $default_comet_params_file .
    """
}

process identification_with_comet {
    maxForks 1
    container 'medbioinf/ident-comparison-comet:latest'

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
workflow comet_identification {
    take:
        sdrf
        fasta
        //mzml
        max_variable_mods

    main:
        default_comet_params_file = create_default_comet_params_file()
        comet_param_files = create_comet_params_file_from_sdrf(sdrf, default_comet_params_file, max_variable_mods)
        //mzidentmls = identification_with_comet(comet_param_file, fasta, mzml)
    emit:
        comet_param_files
}