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
params.mzml_files = '*.mzML'     // may contain globs
params.fasta = 'sprot.fasta'
params.out = './pipeline_results'
params.precursor_tol_ppm = 10
params.fragment_tol_da = 0.02

// default parameter files
params.comet_params_file = "${baseDir}/config/comet.params"
params.maxquant_params_file = "${baseDir}/config/mqpar.xml"
params.msamanda_config_file = "${baseDir}/config/msamanda_settings.xml"
params.msgfplus_params_file = "${baseDir}/config/MSGFPlus_Params.txt"
params.sage_config_file = "${baseDir}/config/sage_config.json"
params.xtandem_config_file = "${baseDir}/config/xtandem_input.xml"

workflow {
    //thermo_raw_files = Channel.fromPath(params.raw_files).flatten()
    mzmls = Channel.fromPath(params.mzml_files).flatten()
    target_fasta = Channel.fromPath(params.fasta).first()

    // get the (default) params / config files
    comet_params_file = Channel.fromPath(params.comet_params_file).first()
    maxquant_params_file = Channel.fromPath(params.maxquant_params_file).first()
    msamanda_config_file = Channel.fromPath(params.msamanda_config_file).first()
    msgfplus_params_file = Channel.fromPath(params.msgfplus_params_file).first()
    sage_config_file = Channel.fromPath(params.sage_config_file).first()
    xtandem_config_file = Channel.fromPath(params.xtandem_config_file).first()

    // create the target-decoy-DB
    target_decoy_fasta = create_decoy_database(target_fasta, "reverse")

    // Identification
    comet_results = comet_identification(comet_params_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    maxquant_results = maxquant_identification(maxquant_params_file, target_fasta, mzmls, params.precursor_tol_ppm)
    msamanda_results = msamanda_identification(msamanda_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    msgfplus_results = msgfplus_identification(msgfplus_params_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm)
    sage_results = sage_identification(sage_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    xtandem_results = xtandem_identification(xtandem_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)

    publish:
    // publish the data in the "out" directory
    comet_results >> 'comet'
    maxquant_results >> 'maxquant'
    msamanda_results >> 'msamanda'
    msgfplus_results >> 'msgfplus'
    sage_results >> 'sage'
    xtandem_results >> 'xtandem'
}

output {
    directory "${params.out}"
    mode 'copy'
}
