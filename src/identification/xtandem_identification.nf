nextflow.enable.dsl=2

/**
 * Creates a X!TAndem params XML file for each RAW file ionthe SDRF file, together with one taxonomy XML file
 * @param sdrf The SDRF file
 * @param sdrf The FASTA file
 * @param max_missed_clavages maximum number of missed cleavages
 * @param max_parent_charge  maximum parent charge

 * @return The XTandem params for each file in the SDRF and the according taxonomy file
 */
process create_xtandem_params_files_from_sdrf {
    container 'medbioinf/ident-comparison-python:latest'

    input:
    path sdrf
    path fasta
    val max_missed_clavages
    val max_parent_charge

    output:
    path "*.tandem_input.xml"
    path "sdrf_convert_taxonomy.xml"

    /** TODO: set the max_threads someway useful */
    """
    python -m sdrf_convert $sdrf tandem --fasta $fasta --max-missed-cleavages $max_missed_clavages --max-parent-charge $max_parent_charge --max-threads 16
    
    # rename absolute paths to current path, to allow for clean passing on in workflow
    workDir=\$(pwd)
    for file in *.tandem_input.xml; do
        sed -i "s;\$workDir/;;g" \$file
    done

    #sed -i "s;\$workDir/;;g" sdrf_convert_taxonomy.xml
    """
}


/**
 * Performs the identifications with XTandem
 */
process identification_with_xtandem {
    maxForks 1
    container 'medbioinf/ident-comparison-xtandem:latest'

    input:
    path xtandem_param_file
    path taxonomy_file
    path fasta
    path mzmls

    output:
    path "*.t.xml"

    """
    tandem $xtandem_param_file
    """
}


/**
 * Exports the identification using Comet configured by a SDRF files
 */
workflow xtandem_identification {
    take:
        sdrf
        fasta
        mzmls
        max_missed_cleavages
        max_parent_charge

    main:
        (xtandem_param_files, taxonomy_file) = create_xtandem_params_files_from_sdrf(sdrf, fasta, max_missed_cleavages, max_parent_charge)
        xtandem_param_files = xtandem_param_files.flatten()

        tandem_xmls = identification_with_xtandem(xtandem_param_files, taxonomy_file, fasta, mzmls.collect())
        tandem_xmls = tandem_xmls.flatten()
    
    emit:
        tandem_xmls
}