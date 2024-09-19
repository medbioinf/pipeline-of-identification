nextflow.enable.dsl=2

xtandem_image = 'quay.io/medbioinf/xtandem:2017.2.1.4'
python_image = 'medbioinf/ident-comparison-python'

// number of threads used by xtandem
params.xtandem_threads = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

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

        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(tandem_xmls, 'xtandem', 'xtandem')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'xtandem')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.drop('xtandem_identification_'.length()).take(it.name.lastIndexOf('.t.xml') - 'xtandem_identification_'.length()) ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'xtandem')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]
    
    emit:
        tandem_xmls
        psm_tsvs
        tda_results
        pout_files
        mokapot_results
        mokapot_features
}


/**
 * Creates a X!TAndem params XML file for each RAW file ionthe SDRF file, together with one taxonomy XML file
 * @param sdrf The SDRF file
 * @param sdrf The FASTA file
 * @param max_missed_clavages maximum number of missed cleavages
 * @param max_parent_charge  maximum parent charge

 * @return The XTandem params for each file in the SDRF and the according taxonomy file
 */
process create_xtandem_params_files_from_sdrf {
    container { python_image }

    input:
    path sdrf
    path fasta
    val max_missed_clavages
    val max_parent_charge

    output:
    path "*.tandem_input.xml"
    path "sdrf_convert_taxonomy.xml"

    """
    python -m sdrf_convert $sdrf tandem --fasta $fasta --max-missed-cleavages $max_missed_clavages --max-parent-charge $max_parent_charge --max-threads $params.xtandem_threads
    
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
    cpus { params.xtandem_threads }
    container { xtandem_image }

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
