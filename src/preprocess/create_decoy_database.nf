nextflow.enable.dsl=2

params.openms_image = 'quay.io/medbioinf/openms:3.1.0'

// number of threads used by maxquant
params.decoy_database_threads = 4


/**
 * Creates a concatenated target-decoy database
 *
 * @return the msms.txt for each mzML file
 */
workflow create_decoy_database {
    take:
    fasta
    decoy_method

    main:
    decoy_fasta = call_decoy_database(fasta, decoy_method)
        
    emit:
    decoy_fasta
}


process call_decoy_database {
    cpus { params.decoy_database_threads }
    container { params.openms_image }

    input:
    path fasta
    val decoy_method

    output:
    path "${fasta.baseName}-rev_decoy.fasta"

    script:
    """
    DecoyDatabase -in ${fasta} -out ${fasta.baseName}-rev_decoy.fasta -decoy_string 'DECOY_' -method ${decoy_method} -threads ${params.decoy_database_threads}
    """
}
