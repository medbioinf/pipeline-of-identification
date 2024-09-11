nextflow.enable.dsl=2

msgfplus_image = 'quay.io/medbioinf/msgfplus:v2024.03.26'

// number of threads used by msamanda
params.msgfplus_threads = 16
params.msgfplus_mem_gb = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'

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

    main:
        msgfplus_results = identification_with_msgfplus(msgfplus_params_file, fasta, mzmls)

        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msgfplus_results, 'mzid', 'msgfplus')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        // percolator
        // ms2rescore
        
    emit:
        msgfplus_results
        psm_tsvs
        pin_files
}


process identification_with_msgfplus {
    cpus { params.msgfplus_threads }
    memory { params.msgfplus_mem_gb + " GB" }
    container { msgfplus_image }

    input:
    path msgfplus_params_file
    path fasta
    path mzmls

    output:
    path "${mzmls.baseName}.mzid"

    """
    java -Xmx${params.msgfplus_mem_gb}G -jar /opt/msgfplus/MSGFPlus.jar -conf ${msgfplus_params_file} -s ${mzmls} -d ${fasta} -thread ${params.msgfplus_threads} -o ${mzmls.baseName}.mzid
    """
}
