nextflow.enable.dsl=2

// parameters set by the command line
//params.raws = 'directory with raw files'
params.mzmls_path = 'directory with mzML files'
params.fasta = 'fasta file'
params.sdrf = 'sdrf file'
params.out = 'output directory'

// TODO: set these in a useful way
params.max_missed_cleavages = 2
params.max_parent_charge = 4

/** TODO: set the num_threads someway useful */

// including modules
include {raw_file_conversion} from workflow.projectDir + '/conversion/raw_file_conversion.nf'
include {xtandem_identification} from workflow.projectDir + '/identification/xtandem_identification.nf'
include {comet_identification} from workflow.projectDir + '/identification/comet_identification.nf'

workflow {
    //thermo_raw_files = Channel.fromPath(params.raws + '/*.raw')
    mzmls = Channel.fromPath(params.mzmls_path + '/*.mzML')
    sdrf = Channel.fromPath(params.sdrf).first()
    fasta = Channel.fromPath(params.fasta).first()

    // Convert raw files to mzML
    //mzmls = raw_file_conversion(thermo_raw_files)

    // Identification
    xtandem_xmls = xtandem_identification(sdrf, fasta, mzmls, params.max_missed_cleavages, params.max_parent_charge)
    comet_tsvs = comet_identification(sdrf, fasta, mzmls)
    
    xtandem_xmls.view { "tandem: $it" }
    comet_tsvs.view { "comet: $it" }
    // Postprocessing

    // Analysis of the data
}