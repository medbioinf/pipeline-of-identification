nextflow.enable.dsl=2

// including modules
include {raw_file_conversion} from workflow.projectDir + '/src/conversion/raw_file_conversion.nf'

include {xtandem_identification} from workflow.projectDir + '/src/identification/xtandem_identification.nf'
include {comet_identification} from workflow.projectDir + '/src/identification/comet_identification.nf'
include {sage_identification} from workflow.projectDir + '/src/identification/sage_identification.nf'

include {pia_tda_analysis} from workflow.projectDir + '/src/postprocessing/pia_tda.nf'


// parameters set by the command line
//params.raws = 'directory with raw files'
params.mzml_files = '*.mzML'     // my contain globs
params.fasta = 'fasta file'
params.sdrf = 'sdrf file'
params.out = 'output directory'

// TODO: set these in a useful way
params.max_missed_cleavages = 2
params.max_parent_charge = 4

/** TODO: set the num_threads someway useful */

// TODO: create these with sdrf-convert
params.sage_config_file = "${baseDir}/config/sage_config.json"

workflow {
    //thermo_raw_files = Channel.fromPath(params.raws + '/*.raw')
    mzmls = Channel.fromPath(params.mzml_files).flatten()
    sdrf = Channel.fromPath(params.sdrf).first()
    fasta = Channel.fromPath(params.fasta).first()

    // TODO: this should go into sdrf-convert
    sage_config_file = Channel.fromPath(params.sage_config_file).first()


    // Convert raw files to mzML
    //mzmls = raw_file_conversion(thermo_raw_files)


    // Identification
    xtandem_xmls = xtandem_identification(sdrf, fasta, mzmls, params.max_missed_cleavages, params.max_parent_charge)
    comet_mzids = comet_identification(sdrf, fasta, mzmls)
    sage_results = sage_identification(sage_config_file, fasta, mzmls)
    
    // Postprocessing

    //pia_tda_comet = pia_tda_analysis(comet_mzids)

    // Analysis of the data
}