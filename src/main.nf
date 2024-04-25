nextflow.enable.dsl=2


params.raws = 'directory with raw fastq files'
params.fasta = 'fasta file'
params.sdrf = 'sdrf file'
params.out = 'output directory'


include {raw_file_conversion} from workflow.projectDir + '/conversion/raw_file_conversion.nf'
include {convert_raw_files_to_mzml} from workflow.projectDir + '/identification/comet_identification.nf'

workflow {
    thermo_raw_files = Channel.fromPath(params.raws + '/*.raw')
    // Convert raw files to mzML
    raw_file_conversion(thermo_raw_files)

    // Identification

    // Preprocessing

    // Analysis of the data

}