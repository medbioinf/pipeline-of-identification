nextflow.enable.dsl=2


params.raws = 'directory with raw fastq files'
params.fasta = 'fasta file'
params.sdrf = 'sdrf file'
params.out = 'output directory'


include {raw_file_conversion} from workflow.projectDir + '/conversion/raw_file_conversion.nf'
include {comet_identification} from workflow.projectDir + '/identification/comet_identification.nf'

workflow {
    thermo_raw_files = Channel.fromPath(params.raws + '/*.raw')
    sdrf = Channel.fromPath(params.sdrf).first()
    fasta = Channel.fromPath(params.fasta).first()
    // Convert raw files to mzML
    mzmls = raw_file_conversion(thermo_raw_files)

    // Identification
    mzidents = comet_identification(sdrf, fasta, mzmls)

    // Preprocessing

    // Analysis of the data

}