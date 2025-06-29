nextflow.enable.dsl=2

params.msamanda_image = 'quay.io/medbioinf/msamanda:3.0.22.071'

// number of threads used by msamanda
params.msamanda_threads = 16
params.msamanda_mem = "64 GB"

params.msamanda_psm_id_pattern = "(.*)"
params.msamanda_spectrum_id_pattern = '(.*)'

include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as onlybest_percolator; psm_percolator as ms2rescore_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'

// msamanda needs explicit "scan=" in the id of a scan (not there in e.g. TimsTOF converted mzML data)
// 1) ms-convert with "--noindex"
// 2) sed -i -e 's;<spectrum\(.*\) id="index=\([0-9]*\)\(.*\);<spectrum\1 id="index=\2 scan=\2\3;' filtered.mzML

/**
 * Executes the identification using MSAmanda
 */
workflow msamanda_identification {
    take:
    msamanda_config_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da

    main:
    msamanda_results = identification_with_msamanda(msamanda_config_file, fasta, mzmls, precursor_tol_ppm, fragment_tol_da)

    // transpose to tuples containing [pin, tsv] files for each mzML 
    return_files = msamanda_results.msamanda_pin.collect()
        .concat(msamanda_results.msamanda_mzid.collect())
        .toList()
        .transpose()
    
    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msamanda_results.msamanda_mzid, 'mzid', 'msamanda')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file
    onlybest_pin_files = psm_tsvs_and_pin.onlybest_pin_file

    pout_files = psm_percolator(pin_files)
    onlybest_pout_files = onlybest_percolator(onlybest_pin_files)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_output.mzid')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.msamanda_psm_id_pattern, params.msamanda_spectrum_id_pattern, '^DECOY_', 'msamanda')

    // perform percolation on MS2Rescore results (both all and onlybest)
    ms2rescore_percolator_results = ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins.concat(ms2rescore_pins.ms2rescore_onlybest_base_pins))

    publish:
    return_files >> 'msamanda'
    psm_tsvs >> 'msamanda'
    pin_files >> 'msamanda'
    onlybest_pin_files >> 'msamanda'
    pout_files >> 'msamanda'
    onlybest_pout_files >> 'msamanda'
    ms2rescore_pins >> 'msamanda'
    ms2rescore_percolator_results >> 'msamanda'
}


process identification_with_msamanda {
    cpus { params.msamanda_threads }
    memory { params.msamanda_mem }
    container { params.msamanda_image }

    input:
    path msamanda_config_file
    path fasta
    path mzmls
    val precursor_tol_ppm
    val fragment_tol_da

    output:
    path "*.mzid_pin.tsv", emit: msamanda_pin
    path "*.mzid", emit: msamanda_mzid

    script:
    """
    cp ${msamanda_config_file} adjusted_msamanda_settings.xml
    sed -i 's;<MS1Tol[^<]*</MS1Tol>;<MS1Tol Unit="ppm">${precursor_tol_ppm}</MS1Tol>;' adjusted_msamanda_settings.xml
    sed -i 's;<MS2Tol[^<]*</MS2Tol>;<MS2Tol Unit="Da">${fragment_tol_da}</MS2Tol>;' adjusted_msamanda_settings.xml

    # MSAmanda command line arguments:
    # Required: -s spectrumFile     single .mgf or .mzml file, or folder with multiple .mgf and .mzml files
    # Required: -d proteinDatabase  single .fasta file or folder with multiple .fasta files, which will be combined into one
    # Required: -e settings.xml
    # Optional: -f fileformat       choose 1 for .csv and 2 for .mzid, default value is 1
    # Optional: -o outputfilename   file or folder where the output should be saved, default path is location of Spectrum file

    MSAmanda -s ${mzmls} -d ${fasta} -e adjusted_msamanda_settings.xml -f 2
    gunzip *.mzid.gz
    """
}
