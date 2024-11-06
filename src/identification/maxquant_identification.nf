nextflow.enable.dsl=2

params.maxquant_image = 'medbioinf/maxquant:v2.6.1.0'

// number of threads used by maxquant
params.maxquant_threads = 16
params.maxquant_mem = "32 GB"

include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from '../postprocessing/default_target_decoy_approach.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'

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
    precursor_tol_ppm

    main:
    maxquant_results = identification_with_maxquant(maxquant_params_file, fasta, mzmls, precursor_tol_ppm)

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(maxquant_results, 'msms', 'maxquant')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    tda_results = target_decoy_approach(psm_tsvs, 'maxquant')

    pout_files = psm_percolator(pin_files)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_msms')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'maxquant')
    ms2rescore_percolator_results = ms2rescore_percolator(ms2rescore_pins)

    publish:
    maxquant_results >> 'maxquant'
    psm_tsvs >> 'maxquant'
    tda_results >> 'maxquant'
    pin_files >> 'maxquant'
    pout_files >> 'maxquant'
    ms2rescore_pins >> 'maxquant'
    ms2rescore_percolator_results >> 'maxquant'
}


process identification_with_maxquant {
    cpus { params.maxquant_threads }
    memory { params.maxquant_mem }
    container { params.maxquant_image }

    stageInMode 'copy'  // MaxQuant respectively Mono does not support symlinks
    
    input:
    path maxquant_params_file
    path fasta
    path mzmls
    val precursor_tol_ppm

    output:
    path "${mzmls.baseName}_msms.txt"

    script:
    """
    # adjust the mqpar.xml file for our current search and path
    cp ${maxquant_params_file} mqpar_adjusted.xml
    sed -i "s;<numThreads>[^<]*</numThreads>;<numThreads>${params.maxquant_threads}</numThreads>;" mqpar_adjusted.xml
    sed -i "s;<fastaFilePath>[^<]*</fastaFilePath>;<fastaFilePath>${fasta}</fastaFilePath>;" mqpar_adjusted.xml
    sed -i "s;<string>CHANGEME_FILE_PATH</string>;<string>${mzmls}</string>;" mqpar_adjusted.xml
    
    sed -i "s;<firstSearchTol>[^<]*</firstSearchTol>;<firstSearchTol>${precursor_tol_ppm}</firstSearchTol>;" mqpar_adjusted.xml
    sed -i "s;<searchTolInPpm>[^<]*</searchTolInPpm>;<searchTolInPpm>True</searchTolInPpm>;" mqpar_adjusted.xml

    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted.xml --changeFolder mqpar_adjusted_new.xml ./ ./

    # execute the identification
    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted_new.xml

    mv combined/txt/msms.txt ${mzmls.baseName}_msms.txt
    """
}
