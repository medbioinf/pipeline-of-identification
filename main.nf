nextflow.enable.dsl=2

nextflow.preview.output = true

// default python image
params.python_image = 'ghcr.io/medbioinf/pipeline-of-identification:latest'
params.oktoberfest_image = 'medbioinf/oktoberfest'

// parameters set by the command line
params.raw_files = ''
params.mzml_files = ''    // may contain globs
params.fasta = ''
params.fasta_target_decoy = ''
params.precursor_tol_ppm = 10
params.fragment_tol_da = 0.02
params.is_timstof = false
params.entrapment_fold = 0
params.use_only_rank1_psms = true

// keep the (converted) mzML files
params.keep_mzmls = true

// should the search engines be executed?
params.execute_comet = true
params.execute_maxquant = true
params.execute_msamanda = true
params.execute_msfragger = true
params.execute_msgfplus = true
params.execute_sage = true
params.execute_xtandem = true

// default parameter files
params.comet_params_file = "${baseDir}/config/comet.params"
params.maxquant_params_file = "${baseDir}/config/mqpar.xml"
params.msamanda_config_file = "${baseDir}/config/msamanda_settings.xml"
params.msfragger_config_file = "${baseDir}/config/closed_fragger.params"
params.msgfplus_params_file = "${baseDir}/config/MSGFPlus_Params.txt"
params.sage_config_file = "${baseDir}/config/sage_config.json"
params.xtandem_config_file = "${baseDir}/config/xtandem_input.xml"

// including modules
include {create_entrapment_database} from "./src/preprocess/create_entrapment_database.nf"
include {create_decoy_database} from "./src/preprocess/create_decoy_database.nf"
include {convert_to_mzml} from "./src/preprocess/convert_to_mzml.nf"

include {comet_identification} from "./src/identification/comet_identification.nf"
include {maxquant_identification} from "./src/identification/maxquant_identification.nf"
include {msamanda_identification} from "./src/identification/msamanda_identification.nf"
include {msfragger_identification} from "./src/identification/msfragger_identification.nf"
include {msgfplus_identification} from "./src/identification/msgfplus_identification.nf"
include {sage_identification} from "./src/identification/sage_identification.nf"
include {xtandem_identification} from "./src/identification/xtandem_identification.nf"

workflow {
    main:

    if (params.is_timstof && !params.raw_files) {
        error("TimsTOF data needs raw files specified!")
    }

    if (!params.raw_files && !params.mzml_files) {
        error("neither raw file/path, nor mzML files given, please provide at least one")
    }

    if (!params.fasta) {
        error("need to state a FASTA file for identification")
    }

    if (!(params.entrapment_fold instanceof Integer)) {
        error("entrapment fold should be an integer")
    }

    if (params.entrapment_fold < 0) {
        error("entrapment fold needs to be >= 0")
    }

    if (params.fasta_target_decoy && params.entrapment_fold > 0) {
        error("entrapment fold and target-decoy FASTA are mutually exclusive, please provide only one")
    }

    if (params.mzml_files) {
        mzmls = Channel.fromPath(params.mzml_files).flatten()
        mzmls_info = mzmls
    } else {
        mzmls_info = "no mzML given"
    }

    if (params.raw_files) {
        raw_files = Channel.fromPath(params.raw_files).flatten()
        raw_files_info = raw_files
    } else {
        raw_files_info = "no raw spectra given"
    }

    fasta_target = Channel.fromPath(params.fasta).first()
    if (params.fasta_target_decoy) {
        fasta_target_decoy = Channel.fromPath(params.fasta_target_decoy).first()
    } else {
        fasta_target_decoy = "no target-decoy FASTA given"
    }

    // color-codes for output: https://github.com/nextflow-io/nf-schema/blob/ecf159f53d45200ef70920c03e75077b5689a386/plugins/nf-schema/src/main/nextflow/validation/Utils.groovy#L222
    log.info "\033[1;32m" + "Pipeline of Identification" + "\033[0;32m" + """\

    raw files:  ${raw_files_info}
    mzML files: ${mzmls_info}
    target FASTA: ${fasta_target}
    target decoy FASTA: ${fasta_target_decoy}
    precursor tolerance (ppm): ${params.precursor_tol_ppm}
    fragment tolerance (Da): ${params.fragment_tol_da}
    timsTOF data: ${params.is_timstof}
    Comet:     ${params.execute_comet}
    MaxQuant:  ${params.execute_maxquant}
    MS Amanda: ${params.execute_msamanda}
    MSFragger: ${params.execute_msfragger}
    MS-GF+:    ${params.execute_msgfplus}
    Sage:      ${params.execute_sage}
    X!Tandem:  ${params.execute_xtandem}
""".stripIndent(true) + "\033[0m"

    // TODO: convert raw files, if not given
    if (!params.mzml_files) {
        mzmls = convert_to_mzml(raw_files)
    }

    if (params.entrapment_fold > 0) {
        // create the entrapment database
        fasta_target = create_entrapment_database(fasta_target, params.entrapment_fold)
    }

    // create the target-decoy-DB, if not given
    if (params.fasta_target_decoy) {
        fasta_target_decoy = Channel.fromPath(params.fasta_target_decoy).first()
    } else {
        fasta_target_decoy = create_decoy_database(fasta_target, "reverse")
    }

    // run the search engines with post-processing
    if (params.execute_comet) {
        comet_params_file = Channel.fromPath(params.comet_params_file).first()
        comet_identification(comet_params_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_maxquant) {
        maxquant_params_file = Channel.fromPath(params.maxquant_params_file).first()
        maxquant_identification(maxquant_params_file, fasta_target, mzmls, params.precursor_tol_ppm)
    }

    if (params.execute_msamanda) {
        msamanda_config_file = Channel.fromPath(params.msamanda_config_file).first()
        msamanda_identification(msamanda_config_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_msfragger) {
        msfragger_config_file = Channel.fromPath(params.msfragger_config_file).first()
        msfragger_identification(msfragger_config_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_msgfplus) {
        msgfplus_params_file = Channel.fromPath(params.msgfplus_params_file).first()
        msgfplus_identification(msgfplus_params_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm)
    }

    if (params.execute_sage) {
        sage_config_file = Channel.fromPath(params.sage_config_file).first()
        sage_identification(sage_config_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }

    if (params.execute_xtandem) {
        xtandem_config_file = Channel.fromPath(params.xtandem_config_file).first()
        xtandem_identification(xtandem_config_file, fasta_target_decoy, mzmls, params.precursor_tol_ppm, params.fragment_tol_da)
    }
}

output {
    'mzmls' {
        enabled params.keep_mzmls
        path 'mzmls'
    }

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

    'msfragger' {
        enabled params.execute_msfragger
        path 'msfragger'
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