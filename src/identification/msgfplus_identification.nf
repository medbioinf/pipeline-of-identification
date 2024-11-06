nextflow.enable.dsl=2

params.msgfplus_image = 'quay.io/medbioinf/msgfplus:v2024.03.26'

// number of threads used by msamanda
params.msgfplus_threads = 16
params.msgfplus_mem_gb = 16

include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from '../postprocessing/default_target_decoy_approach.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'

/**
 * Executes the identification using MS-GF+
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow msgfplus_identification {
    take:
    msgfplus_params_file
    fasta
    mzmls
    precursor_tol_ppm

    main:
    msgfplus_results = identification_with_msgfplus(msgfplus_params_file, fasta, mzmls, precursor_tol_ppm)

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msgfplus_results, 'mzid', 'msgfplus')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    tda_results = target_decoy_approach(psm_tsvs, 'msgfplus')

    pout_files = psm_percolator(pin_files)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'msgfplus')
    ms2rescore_percolator_results = ms2rescore_percolator(ms2rescore_pins)
        
    publish:
    msgfplus_results >> 'msgfplus'
    psm_tsvs >> 'msgfplus'
    tda_results >> 'msgfplus'
    pin_files >> 'msgfplus'
    pout_files >> 'msgfplus'
    ms2rescore_pins >> 'msgfplus'
    ms2rescore_percolator_results >> 'msgfplus'
}


process identification_with_msgfplus {
    cpus { params.msgfplus_threads }
    memory { params.msgfplus_mem_gb + " GB" }
    container { params.msgfplus_image }

    input:
    path msgfplus_params_file
    path fasta
    path mzmls
    val precursor_tol_ppm

    output:
    path "${mzmls.baseName}.mzid"

    script:
    """
    cp ${msgfplus_params_file} adjusted_MSGFPlus_Params.txt
    sed -i 's;PrecursorMassTolerance=.*;PrecursorMassTolerance=${precursor_tol_ppm};' adjusted_MSGFPlus_Params.txt

    java -Xmx${params.msgfplus_mem_gb}G -jar /opt/msgfplus/MSGFPlus.jar -conf adjusted_MSGFPlus_Params.txt -s ${mzmls} -d ${fasta} -thread ${params.msgfplus_threads} -o ${mzmls.baseName}.mzid
    """
}
