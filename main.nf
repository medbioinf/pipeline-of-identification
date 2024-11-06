nextflow.enable.dsl=2

nextflow.preview.output = true

// default python image
params.python_image = 'medbioinf/ident-comparison-python'

// parameters set by the command line
params.mzml_files = '*.mzML'    // may contain globs
params.fasta = 'sprot.fasta'    // single file
params.precursor_tol_ppm = 10
params.fragment_tol_da = 0.02

// should the search engines be executed?
params.execute_comet = true
params.execute_maxquant = true
params.execute_msamanda = true
params.execute_msgfplus = true
params.execute_sage = true
params.execute_xtandem = true

// default parameter files
params.comet_params_file = "${baseDir}/config/comet.params"
params.maxquant_params_file = "${baseDir}/config/mqpar.xml"
params.msamanda_config_file = "${baseDir}/config/msamanda_settings.xml"
params.msgfplus_params_file = "${baseDir}/config/MSGFPlus_Params.txt"
params.sage_config_file = "${baseDir}/config/sage_config.json"
params.xtandem_config_file = "${baseDir}/config/xtandem_input.xml"

// including modules
include {create_decoy_database} from "./src/preprocess/create_decoy_database.nf"

include {comet_identification} from "./src/identification/comet_identification.nf"
include {maxquant_identification} from "./src/identification/maxquant_identification.nf"
include {msamanda_identification} from "./src/identification/msamanda_identification.nf"
include {msgfplus_identification} from "./src/identification/msgfplus_identification.nf"
include {sage_identification} from "./src/identification/sage_identification.nf"
include {xtandem_identification} from "./src/identification/xtandem_identification.nf"


workflow {
    main:
    mzmls = Channel.fromPath(params.mzml_files).flatten()
    target_fasta = Channel.fromPath(params.fasta).first()

    // create the target-decoy-DB
    target_decoy_fasta = create_decoy_database(target_fasta, "reverse")

    // TODO: checken, welche params sich in "high res / low res" unterscheiden
    // - auf jeden fall comet
    // - msgf+ auch das modell

    // TODO: add PTM definitions

    // run the search engines with post-processing
    if (params.execute_comet) {
        comet_params_file = Channel.fromPath(params.comet_params_file).first()
        comet_identification(comet_params_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_maxquant) {
        maxquant_params_file = Channel.fromPath(params.maxquant_params_file).first()
        maxquant_identification(maxquant_params_file, target_fasta, mzmls, params.precursor_tol_ppm)
    }

    if (params.execute_msamanda) {
        msamanda_config_file = Channel.fromPath(params.msamanda_config_file).first()
        msamanda_identification(msamanda_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_msgfplus) {
        msgfplus_params_file = Channel.fromPath(params.msgfplus_params_file).first()
        msgfplus_identification(msgfplus_params_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm)
    }

    if (params.execute_sage) {
        sage_config_file = Channel.fromPath(params.sage_config_file).first()
        sage_identification(sage_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_xtandem) {
        xtandem_config_file = Channel.fromPath(params.xtandem_config_file).first()
        xtandem_identification(xtandem_config_file, target_decoy_fasta, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }
}

output {
    'comet' {
        enabled params.execute_comet
        path 'comet'
    }

    'maxquant' {
        enabled params.execute_maxquant
        path 'maxquant'
    }

    'msamanda' {
        enabled params.execute_msamanda
        path 'msamanda'
    }

    'msgfplus' {
        enabled params.execute_msgfplus
        path 'msgfplus'
    }

    'sage' {
        enabled params.execute_sage
        path 'sage'
    }

    'xtandem' {
        enabled params.execute_xtandem
        path 'xtandem'
    }
}