nextflow.enable.dsl=2

msamanda_image = 'quay.io/medbioinf/msamanda:3.0.22.071'

// number of threads used by msamanda
params.msamanda_threads = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

/**
 * Executes the identification using Sage
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow msamanda_identification {
    take:
        msamanda_config_file
        fasta
        mzmls

    main:
        msamanda_results = identification_with_msamanda(msamanda_config_file, fasta, mzmls)

        // transpose to tuples containing [pin, tsv] files for each mzML 
        return_files = msamanda_results.msamanda_pin.collect()
            .concat(msamanda_results.msamanda_mzid.collect())
            .toList()
            .transpose()
        
        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msamanda_results.msamanda_mzid, 'mzid', 'msamanda')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'msamanda')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_output.mzid')) + '.mzML'  ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'msamanda')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]

    emit:
        return_files
        psm_tsvs
        tda_results
        pin_files
        pout_files
        mokapot_results
        mokapot_features
}


process identification_with_msamanda {
    cpus { params.msamanda_threads }
    container { msamanda_image }

    input:
    path msamanda_config_file
    path fasta
    path mzmls

    output:
    path "*.mzid_pin.tsv", emit: msamanda_pin
    path "*.mzid", emit: msamanda_mzid

    """
    # Required: -s spectrumFile     single .mgf or .mzml file, or folder with multiple .mgf and .mzml files
    # Required: -d proteinDatabase  single .fasta file or folder with multiple .fasta files, which will be combined into one
    # Required: -e settings.xml
    # Optional: -f fileformat       choose 1 for .csv and 2 for .mzid, default value is 1
    # Optional: -o outputfilename   file or folder where the output should be saved, default path is location of Spectrum file

    MSAmanda -s ${mzmls} -d ${fasta} -e ${msamanda_config_file} -f 2
    gunzip *.mzid.gz
    """
}
