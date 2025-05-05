nextflow.enable.dsl=2

params.msgfplus_image = 'quay.io/medbioinf/msgfplus:v2024.03.26'

// params for MS-GF+
params.msgfplus_threads = 6
params.msgfplus_mem_gb = 16
params.msgfplus_tasks = 0

params.msgfplus_instrument = "1" // 0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR/Lumos, 2: TOF, 3: Q-Exactive

params.msgfplus_split_input = 10000     // split input mzMLs into chunks of this size, 0 to disable

params.msgfplus_psm_id_pattern = "(.*)"
params.msgfplus_spectrum_id_pattern = '(.*)'

include {convert_chunked_result_to_psm_utils; enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as onlybest_percolator; psm_percolator as ms2rescore_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'

include {split_mzml_into_chunks} from '../preprocess/convert_to_mzml.nf'

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
    fasta_index = build_msgfplus_index(fasta)

    if (params.msgfplus_split_input > 0) {
        chunked_mzmls = split_mzml_into_chunks(params.msgfplus_split_input, mzmls)
        mzmls_to_chunks = chunked_mzmls.transpose()
    } else {
        mzmls_to_chunks = mzmls.map{ it -> [it.baseName, it] }
    }

    msgfplus_results = identification_with_msgfplus(msgfplus_params_file, fasta_index, mzmls_to_chunks, precursor_tol_ppm)

    psm_tsvs_with_mzml = convert_chunked_result_to_psm_utils(msgfplus_results, 'mzid')

    grouped_results = psm_tsvs_with_mzml.groupTuple(by: 0)
    if (params.msgfplus_split_input > 0) {
        merged_results = merge_psms(grouped_results)
    } else {
        merged_results = psm_tsvs_with_mzml.map{ it -> it[1] }
    }

    psm_tsvs_and_pin = enhance_psm_tsv(merged_results, 'msgfplus')

    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file
    onlybest_pin_files = psm_tsvs_and_pin.onlybest_pin_file

    pout_files = psm_percolator(pin_files)
    onlybest_pout_files = onlybest_percolator(onlybest_pin_files)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.msgfplus_psm_id_pattern, params.msgfplus_spectrum_id_pattern, '^DECOY_', 'msgfplus')

    // perform percolation on MS2Rescore results (both all and onlybest)
    ms2rescore_percolator_results = ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins.concat(ms2rescore_pins.ms2rescore_onlybest_base_pins))

    publish:
    msgfplus_results.map{ it -> it[1] } >> 'msgfplus'
    psm_tsvs >> 'msgfplus'
    pin_files >> 'msgfplus'
    onlybest_pin_files >> 'msgfplus'
    pout_files >> 'msgfplus'
    onlybest_pout_files >> 'msgfplus'
    ms2rescore_pins >> 'msgfplus'
    ms2rescore_percolator_results >> 'msgfplus'
}

process identification_with_msgfplus {
    cpus { params.msgfplus_threads }
    memory { params.msgfplus_mem_gb + " GB" }
    container { params.msgfplus_image }

    input:
    path msgfplus_params_file
    tuple path(fasta), path(canno), path(cnlcp), path(csarr), path(cseq)
    tuple val(original_mzml_basename), path(mzml)
    val precursor_tol_ppm

    output:
    tuple val(original_mzml_basename), path("${mzml.baseName}.mzid")

    script:
    """
    cp ${msgfplus_params_file} adjusted_MSGFPlus_Params.txt
    sed -i 's;^PrecursorMassTolerance=.*;PrecursorMassTolerance=${precursor_tol_ppm};' adjusted_MSGFPlus_Params.txt
    sed -i 's;^InstrumentID=.*;InstrumentID=${params.msgfplus_instrument};' adjusted_MSGFPlus_Params.txt

    java -Xmx${params.msgfplus_mem_gb}G -jar /opt/msgfplus/MSGFPlus.jar -conf adjusted_MSGFPlus_Params.txt -s ${mzml} -d ${fasta} -thread ${params.msgfplus_threads} -tasks ${params.msgfplus_tasks} -o ${mzml.baseName}.mzid
    """
}


process build_msgfplus_index {
    cpus { params.msgfplus_threads }
    memory { params.msgfplus_mem_gb + " GB" }
    container { params.msgfplus_image }

    input:
    path fasta

    output:
    tuple path("${fasta}"), path("${fasta.baseName}.canno"), path("${fasta.baseName}.cnlcp"), path("${fasta.baseName}.csarr"), path("${fasta.baseName}.cseq")

    script:
    """
    java -Xmx${params.msgfplus_mem_gb}G -cp /opt/msgfplus/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d ${fasta} -tda 0 -o ./ -decoy DECOY_
    """
}


process merge_psms {
    cpus 2
    memory { params.msgfplus_mem_gb + " GB" }
    container { params.python_image }

    input:
    tuple val(original_mzml_basename), path(psm_tsvs)

    output:
    path "${original_mzml_basename}.mzid.psm_utils.tsv"

    script:
    """
    merge_chunked_psm_files.py --org_filebase ${original_mzml_basename} --out_filename ${original_mzml_basename}.mzid.psm_utils.tsv --files ${psm_tsvs}
    """
}