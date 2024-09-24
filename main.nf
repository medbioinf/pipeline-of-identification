nextflow.enable.dsl=2

nextflow.preview.output = true

// including modules
include {raw_file_conversion} from workflow.projectDir + '/src/preprocess/raw_file_conversion.nf'
include {create_decoy_database} from workflow.projectDir + '/src/preprocess/create_decoy_database.nf'

include {xtandem_identification} from workflow.projectDir + '/src/identification/xtandem_identification.nf'
include {comet_identification} from workflow.projectDir + '/src/identification/comet_identification.nf'
include {sage_identification} from workflow.projectDir + '/src/identification/sage_identification.nf'
include {msamanda_identification} from workflow.projectDir + '/src/identification/msamanda_identification.nf'
include {msgfplus_identification} from workflow.projectDir + '/src/identification/msgfplus_identification.nf'
include {maxquant_identification} from workflow.projectDir + '/src/identification/maxquant_identification.nf'


// parameters set by the command line
//params.raws = 'directory with raw files'
params.mzml_files = '*.mzML'     // my contain globs
params.fasta = 'fasta file (without decoys!)'
params.sdrf = 'sdrf file'
params.out = 'output directory'

// TODO: set these in a useful way
params.max_missed_cleavages = 2
params.max_parent_charge = 4

/** TODO: set the num_threads someway useful */

// TODO: create these with sdrf-convert
params.sage_config_file = "${baseDir}/config/sage_config.json"
params.msamanda_config_file = "${baseDir}/config/msamanda_settings.xml"
params.msgfplus_params_file = "${baseDir}/config/MSGFPlus_Params.txt"
params.maxquant_params_file = "${baseDir}/config/mqpar.xml"

workflow {
    //thermo_raw_files = Channel.fromPath(params.raw_files).flatten()
    mzmls = Channel.fromPath(params.mzml_files).flatten()
    sdrf = Channel.fromPath(params.sdrf).first()
    target_fasta = Channel.fromPath(params.fasta).first()

    // TODO: this should go into sdrf-convert
    sage_config_file = Channel.fromPath(params.sage_config_file).first()
    msamanda_config_file = Channel.fromPath(params.msamanda_config_file).first()
    msgfplus_params_file = Channel.fromPath(params.msgfplus_params_file).first()
    maxquant_params_file = Channel.fromPath(params.maxquant_params_file).first()

    // Convert raw files to mzML
    //mzmls = raw_file_conversion(thermo_raw_files)

    // create the target-decoy-DB
    target_decoy_fasta = create_decoy_database(target_fasta, "reverse")

    // Identification
    xtandem_results = xtandem_identification(sdrf, target_decoy_fasta, mzmls, params.max_missed_cleavages, params.max_parent_charge)
    comet_mzids = comet_identification(sdrf, target_decoy_fasta, mzmls)
    sage_results = sage_identification(sage_config_file, target_decoy_fasta, mzmls)
    msamanda_results = msamanda_identification(msamanda_config_file, target_decoy_fasta, mzmls)
    msgfplus_results = msgfplus_identification(msgfplus_params_file, target_decoy_fasta, mzmls)
    maxquant_results = maxquant_identification(maxquant_params_file, target_fasta, mzmls)
    
    // Postprocessing

    //pia_tda_comet = pia_tda_analysis(comet_mzids)

    // Analysis of the data

    publish:
    // publish the data in the "out" directory
    xtandem_results >> 'xtandem'
    comet_mzids >> 'comet'
    sage_results >> 'sage'
    msamanda_results >> 'msamanda'
    msgfplus_results >> 'msgfplus'
    maxquant_results >> 'maxquant'
}

output {
    directory "${params.out}"
    mode 'copy'
}
