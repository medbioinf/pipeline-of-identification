nextflow.enable.dsl=2

maxquant_image = 'medbioinf/maxquant'

// number of threads used by maxquant
params.maxquant_threads = 16

/**
 * Executes the identification using MaxQuant
 *
 * @return the msms.txt for each mzML file
 */
workflow maxquant_identification {
    take:
        maxquant_params_file
        fasta
        mzmls

    main:
        maxquant_results = identification_with_maxquant(maxquant_params_file, fasta, mzmls)
        
    emit:
        maxquant_results
}


process identification_with_maxquant {
    cpus { params.maxquant_threads }
    container { maxquant_image }

    input:
    path maxquant_params_file
    path fasta
    path mzmls

    output:
    path "${mzmls.baseName}-msms.txt"

    """
    # adjust the mqpar.xml file for our current search and path
    cp ${maxquant_params_file} mqpar_adjusted.xml
    sed -i "s;<numThreads>[^<]*</numThreads>;<numThreads>${params.maxquant_threads}</numThreads>;" mqpar_adjusted.xml
    sed -i "s;<fastaFilePath>[^<]*</fastaFilePath>;<fastaFilePath>${fasta}</fastaFilePath>;" mqpar_adjusted.xml
    sed -i "s;<string>CHANGEME_FILE_PATH</string>;<string>${mzmls}</string>;" mqpar_adjusted.xml

    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted.xml --changeFolder mqpar_adjusted_new.xml ./ ./

    # execute the identification
    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted_new.xml

    # rename the file to contain the mzML basename
    mv ./combined/txt/msms.txt ./${mzmls.baseName}-msms.txt
    """
}
