nextflow.enable.dsl=2

/**
 * Converts thermo raw files to mzML using msconvert
 */
process convert_thermo_raw_files {
    container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses'

    input:
    path raw

    output:
    path "${raw.baseName}.mzML"


    script:
    """
    wine msconvert ${raw} --mzML --zlib --filter "peakPicking true 1-"
    """
}

/**
 * Exports raw file conversion
 */
workflow raw_file_conversion {
    take:
        thermo_raw_files // a list of thermo raw files
    main:
        thermo_mzmls = convert_thermo_raw_files(thermo_raw_files)
        mzmls = thermo_mzmls
    emit:
        mzmls
}