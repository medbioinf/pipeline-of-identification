nextflow.enable.dsl=2

fdrbench_image = 'quay.io/medbioinf/fdrbench-nightly:146f77'


/**
 * Adds decoys and/or entapments to the FASTA file.
 *
 * @param fasta Path to the FASTA file
 *
 * @return Path to the new FASTA file
 */
workflow create_entrapment_database {
    take:
        fasta
        fold

    main:
        entrapment_fasta = call_entrapment_database(fasta, fold)
        
    emit:
        entrapment_fasta
}


/**
 * Adds entrapments to the FASTA file using FDRBench.
 * https://doi.org/10.1101/2024.06.01.596967
 *
 * @param fasta Path to the FASTA file
 * @param fold Fold change for entrapment
 *
 * @return Path to the new FASTA file
 */
process call_entrapment_database {
    cpus 1
    container { fdrbench_image }

    input: 
    path fasta
    val fold

    output:
    path "${fasta.baseName}-entrapment.fasta"

    script:
    """
    java -jar /opt/fdrbench/fdrbench.jar -db ${fasta} -o ${fasta.baseName}-entrapment.fasta -fold ${fold} -level protein -entrapment_label ENTRAPMENT_ -entrapment_pos 0 -uniprot -check
    # 'Reaheader' to add entrapment index to database and accesion part of the header
    sed -r -i "s;^(>ENTRAPMENT_tr)(\\|.+)(\\|.+\\_(.+))\$;\\1_\\4\\2_ENTR_\\4\\3;g" ${fasta.baseName}-entrapment.fasta
    """

}
