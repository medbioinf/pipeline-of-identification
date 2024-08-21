nextflow.enable.dsl=2

msamanda_image = 'quay.io/medbioinf/msamanda:3.0.22.071'
python_image = 'medbioinf/ident-comparison-python'

// number of threads used by msamanda
params.msamanda_threads = 16


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
        
    emit:
        return_files
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