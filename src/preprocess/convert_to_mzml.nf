nextflow.enable.dsl=2

params.msconvert_image = 'proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25073-842baef'


workflow convert_to_mzml {
    take:
    input_path

    main:
    if (params.is_timstof) {
        mzml = convert_bruker_d(input_path)
    } else {
        mzml = convert_thermo_raw(input_path)
    }
        
    emit:
    mzml

    publish:
    mzml >> 'mzmls'
}

process convert_thermo_raw {
    cpus 2
    container { params.msconvert_image }

    input:
    path input_raw

    output:
    path "${input_raw.baseName}.mzML"

    script:
    """
    wine msconvert --mzML --zlib=off --filter "peakPicking vendor msLevel=1-" ${input_raw} --outfile ${input_raw.baseName}.mzML
    """
}

process convert_bruker_d {
    cpus 2
    container { params.python_image }

    input:
    path input_raw

    output:
    path "${input_raw.baseName}.mzML"

    script:
    """
    # here goes the converter...

    # msamanda needs explicit "scan=" in the id of a scan (not there in e.g. TimsTOF converted mzML data)
    sed -i -e 's;<spectrum\\(.*\\) id="index=\\([0-9]*\\)\\(.*\\);<spectrum\\1 id="index=\\2 scan=\\2\\3;' ${input_raw.baseName}.mzML
    """
}


process split_mzml_into_chunks {
    cpus 2
    memory "8 GB"
    container { params.msconvert_image }

    input:
    val chunksize
    path mzml_file

    output:
    tuple val(mzml_file.baseName), path("${mzml_file.baseName}--*.mzML")

    script:
    """
    CHUNKSIZE=${chunksize}

    # get only MS2 spectra
    wine msconvert --mzML --filter "msLevel 2-" "${mzml_file}" --outfile "${mzml_file.baseName}-MS2only.mzML"

    # get the last spectrum 
    MAX_INDEX=\$(tac ${mzml_file.baseName}-MS2only.mzML | grep -m1 "<spectrum.*index=" | sed 's;.*index="\\([0-9]*\\)".*;\\1;')
    echo "max index: \${MAX_INDEX}"

    for ((i=0; i<=MAX_INDEX; i+=CHUNKSIZE)); do
        START_INDEX=\$i
        END_INDEX=\$((i+CHUNKSIZE-1))
        if ((END_INDEX > MAX_INDEX)); then
            END_INDEX=\$MAX_INDEX
        fi
        echo "Processing spectra from index \${START_INDEX} to \${END_INDEX}"

        wine msconvert --mzML --filter "index [\${START_INDEX},\${END_INDEX}]" "${mzml_file.baseName}-MS2only.mzML" --outfile "${mzml_file.baseName}--\${START_INDEX}_\${END_INDEX}.mzML"
    done
    """
}
