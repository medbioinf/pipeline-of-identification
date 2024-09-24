nextflow.enable.dsl=2

maxquant_image = 'medbioinf/maxquant:v2.6.1.0'

// number of threads used by maxquant
params.maxquant_threads = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

/**
 * Executes the identification using MaxQuant
 *
 * @return the msms.txt for each mzML file
 */
workflow maxquant_identification {
    take:
        maxquant_params_file
        fasta
        mzmls

    main:
        maxquant_results = identification_with_maxquant(maxquant_params_file, fasta, mzmls)

        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(maxquant_results, 'msms', 'maxquant')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'maxquant')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_msms')) + '.mzML'  ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'maxquant')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]

    emit:
        maxquant_results
        psm_tsvs
        tda_results
        pin_files
        pout_files
        mokapot_results
        mokapot_features
}


process identification_with_maxquant {
    cpus { params.maxquant_threads }
    container { maxquant_image }

    stageInMode 'copy'  // MaxQuant respectively Mono does not support symlinks
    
    input:
    path maxquant_params_file
    path fasta
    path mzmls

    output:
    path "${mzmls.baseName}_msms.txt"

    """
    # adjust the mqpar.xml file for our current search and path
    cp ${maxquant_params_file} mqpar_adjusted.xml
    sed -i "s;<numThreads>[^<]*</numThreads>;<numThreads>${params.maxquant_threads}</numThreads>;" mqpar_adjusted.xml
    sed -i "s;<fastaFilePath>[^<]*</fastaFilePath>;<fastaFilePath>${fasta}</fastaFilePath>;" mqpar_adjusted.xml
    sed -i "s;<string>CHANGEME_FILE_PATH</string>;<string>${mzmls}</string>;" mqpar_adjusted.xml

    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted.xml --changeFolder mqpar_adjusted_new.xml ./ ./

    # execute the identification
    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted_new.xml

    mv combined/txt/msms.txt ${mzmls.baseName}_msms.txt
    """
}
