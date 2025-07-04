nextflow.enable.dsl=2

params.maxquant_image = 'medbioinf/maxquant:v2.6.1.0'

// number of threads used by maxquant
params.maxquant_threads = 4
params.maxquant_mem = "32 GB"

params.maxquant_psm_id_pattern = ""
params.maxquant_spectrum_id_pattern = ""

include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'

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
    precursor_tol_ppm

    main:
    // for TimsTOF data, always process the .d path instead of the mzML files
    if (params.is_timstof) {
        process_files = Channel.fromPath(params.raw_files).flatten()
    } else {
        process_files = mzmls
    }

    maxquant_results = identification_with_maxquant(maxquant_params_file, fasta, process_files, precursor_tol_ppm)

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(maxquant_results, 'msms', 'maxquant')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    pout_files = psm_percolator(pin_files)

    if (params.maxquant_psm_id_pattern) {
        psm_id_pattern = params.maxquant_psm_id_pattern
    } else {
        psm_id_pattern = "(.*)"
    }
    if (params.maxquant_spectrum_id_pattern) {
        spectrum_id_pattern = params.maxquant_spectrum_id_pattern
    } else{
        if (params.is_timstof) {
            spectrum_id_pattern = '(.*)'
        } else {
            spectrum_id_pattern = '.*scan=(\\d+)$'
        }
    }

    if (params.is_timstof) {
        psm_tsvs_and_spectrafiles = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_msms')) + '.d'  ] }
    } else {
        psm_tsvs_and_spectrafiles = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_msms')) + '.mzML'  ] }
    }

    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_spectrafiles, psm_tsvs.collect(), process_files.collect(), psm_id_pattern, spectrum_id_pattern, '', 'maxquant')
    
    // perform percolation on MS2Rescore results
    ms2rescore_percolator_results = ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins)

    publish:
    maxquant_results >> 'maxquant'
    psm_tsvs >> 'maxquant'
    pin_files >> 'maxquant'
    pout_files >> 'maxquant'
    ms2rescore_pins >> 'maxquant'
    ms2rescore_percolator_results >> 'maxquant'
}


process identification_with_maxquant {
    cpus { params.maxquant_threads }
    memory { params.maxquant_mem }
    container { params.maxquant_image }

    stageInMode 'copy'  // MaxQuant respectively Mono does not support symlinks
    
    input:
    path maxquant_params_file
    path fasta
    path mzmls
    val precursor_tol_ppm

    output:
    path "${mzmls.baseName}_msms.txt"

    script:
    if (params.is_timstof) {
        maxquant_quantmode = 2
        maxquant_msinstrument = 4
        maxquant_usems1centroids = 'False'
        maxquant_usems2centroids = 'False'
        maxquant_intensitydetermination = 1
        maxquant_advancedpeaksplitting = 'True'
        maxquant_intensitythresholds1dda = 35
        maxquant_lcmsruntype = 'TIMS-DDA'
        maxquant_lfqmode = 1
        maxquant_mainsearchtol = 10
        maxquant_isotopematchtol = 0.005
        maxquant_isotopematchtolinppm = 'False'
        maxquant_checkmassdeficit = 'False'
        maxquant_intensitydependentcalibration = 'True'
        maxquant_minscoreforcalibration = 40
        maxquant_timshalfwidth = 10
        maxquant_timsstep = 3
        maxquant_timsresolution = 35000
        maxquant_timsminmsmsintensity = 1.5
        maxquant_lfqtopncorrelatingpeptides = 100
        maxquant_lfqpeptidecorrelation = 0
    } else {
        maxquant_quantmode = 1
        maxquant_msinstrument = 0
        maxquant_usems1centroids = 'True'
        maxquant_usems2centroids = 'True'
        maxquant_intensitydetermination = 0
        maxquant_advancedpeaksplitting = 'False'
        maxquant_intensitythresholds1dda = 0
        maxquant_lcmsruntype = 'Standard'
        maxquant_lfqmode = 0
        maxquant_mainsearchtol = 4.5
        maxquant_isotopematchtol = 2
        maxquant_isotopematchtolinppm = 'True'
        maxquant_checkmassdeficit = 'True'
        maxquant_intensitydependentcalibration = 'False'
        maxquant_minscoreforcalibration = 70
        maxquant_timshalfwidth = 0
        maxquant_timsstep = 0
        maxquant_timsresolution = 0
        maxquant_timsminmsmsintensity = 0
        maxquant_lfqtopncorrelatingpeptides = 3
        maxquant_lfqpeptidecorrelation = 3
    }
    """
    # adjust the mqpar.xml file for our current search and path
    cp ${maxquant_params_file} mqpar_adjusted.xml

    sed -i "s;<fastaFilePath>[^<]*</fastaFilePath>;<fastaFilePath>${fasta}</fastaFilePath>;" mqpar_adjusted.xml

    sed -i "s;<string>CHANGEME_FILE_PATH</string>;<string>${mzmls}</string>;" mqpar_adjusted.xml

    sed -i "s;<quantMode>[^<]*</quantMode>;<quantMode>${maxquant_quantmode}</quantMode>;" mqpar_adjusted.xml

    sed -i "s;<numThreads>[^<]*</numThreads>;<numThreads>${params.maxquant_threads}</numThreads>;" mqpar_adjusted.xml

    sed -i "s;<msInstrument>[^<]*</msInstrument>;<msInstrument>${maxquant_msinstrument}</msInstrument>;" mqpar_adjusted.xml
    sed -i "s;<useMs1Centroids>[^<]*</useMs1Centroids>;<useMs1Centroids>${maxquant_usems1centroids}</useMs1Centroids>;" mqpar_adjusted.xml
    sed -i "s;<useMs2Centroids>[^<]*</useMs2Centroids>;<useMs2Centroids>${maxquant_usems2centroids}</useMs2Centroids>;" mqpar_adjusted.xml
    sed -i "s;<intensityDetermination>[^<]*</intensityDetermination>;<intensityDetermination>${maxquant_intensitydetermination}</intensityDetermination>;" mqpar_adjusted.xml
    sed -i "s;<advancedPeakSplitting>[^<]*</advancedPeakSplitting>;<advancedPeakSplitting>${maxquant_advancedpeaksplitting}</advancedPeakSplitting>;" mqpar_adjusted.xml
    sed -i "s;<intensityThresholdMs1Dda>[^<]*</intensityThresholdMs1Dda>;<intensityThresholdMs1Dda>${maxquant_intensitythresholds1dda}</intensityThresholdMs1Dda>;" mqpar_adjusted.xml
    sed -i "s;<lcmsRunType>[^<]*</lcmsRunType>;<lcmsRunType>${maxquant_lcmsruntype}</lcmsRunType>;" mqpar_adjusted.xml
    sed -i "s;<lfqMode>[^<]*</lfqMode>;<lfqMode>${maxquant_lfqmode}</lfqMode>;" mqpar_adjusted.xml

    sed -i "s;<firstSearchTol>[^<]*</firstSearchTol>;<firstSearchTol>${precursor_tol_ppm}</firstSearchTol>;" mqpar_adjusted.xml
    sed -i "s;<mainSearchTol>[^<]*</mainSearchTol>;<mainSearchTol>${maxquant_mainsearchtol}</mainSearchTol>;" mqpar_adjusted.xml
    sed -i "s;<searchTolInPpm>[^<]*</searchTolInPpm>;<searchTolInPpm>True</searchTolInPpm>;" mqpar_adjusted.xml
    sed -i "s;<isotopeMatchTol>[^<]*</isotopeMatchTol>;<isotopeMatchTol>${maxquant_isotopematchtol}</isotopeMatchTol>;" mqpar_adjusted.xml
    sed -i "s;<isotopeMatchTolInPpm>[^<]*</isotopeMatchTolInPpm>;<isotopeMatchTolInPpm>${maxquant_isotopematchtolinppm}</isotopeMatchTolInPpm>;" mqpar_adjusted.xml
    sed -i "s;<checkMassDeficit>[^<]*</checkMassDeficit>;<checkMassDeficit>${maxquant_checkmassdeficit}</checkMassDeficit>;" mqpar_adjusted.xml
    sed -i "s;<intensityDependentCalibration>[^<]*</intensityDependentCalibration>;<intensityDependentCalibration>${maxquant_intensitydependentcalibration}</intensityDependentCalibration>;" mqpar_adjusted.xml
    sed -i "s;<minScoreForCalibration>[^<]*</minScoreForCalibration>;<minScoreForCalibration>${maxquant_minscoreforcalibration}</minScoreForCalibration>;" mqpar_adjusted.xml
    sed -i "s;<timsHalfWidth>[^<]*</timsHalfWidth>;<timsHalfWidth>${maxquant_timshalfwidth}</timsHalfWidth>;" mqpar_adjusted.xml
    sed -i "s;<timsStep>[^<]*</timsStep>;<timsStep>${maxquant_timsstep}</timsStep>;" mqpar_adjusted.xml
    sed -i "s;<timsResolution>[^<]*</timsResolution>;<timsResolution>${maxquant_timsresolution}</timsResolution>;" mqpar_adjusted.xml
    sed -i "s;<timsMinMsmsIntensity>[^<]*</timsMinMsmsIntensity>;<timsMinMsmsIntensity>${maxquant_timsminmsmsintensity}</timsMinMsmsIntensity>;" mqpar_adjusted.xml

    sed -i "s;<lfqTopNCorrelatingPeptides>[^<]*</lfqTopNCorrelatingPeptides>;<lfqTopNCorrelatingPeptides>${maxquant_lfqtopncorrelatingpeptides}</lfqTopNCorrelatingPeptides>;" mqpar_adjusted.xml
    sed -i "s;<lfqPeptideCorrelation>[^<]*</lfqPeptideCorrelation>;<lfqPeptideCorrelation>${maxquant_lfqpeptidecorrelation}</lfqPeptideCorrelation>;" mqpar_adjusted.xml

    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted.xml --changeFolder mqpar_adjusted_new.xml ./ ./

    # execute the identification
    dotnet /opt/MaxQuant/bin/MaxQuantCmd.dll mqpar_adjusted_new.xml

    mv combined/txt/msms.txt ${mzmls.baseName}_msms.txt
    """
}
