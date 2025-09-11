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
        entrapment_fasta = call_entrapment_database(fasta, fold, params.fdrbench_mem_gb)
        
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
    memory "${ memory_limit }.GB"

    label 'fdrbench_image'

    input: 
    path fasta
    val fold
    val memory_limit

    output:
    path "${fasta.baseName}-entrapment.fasta"

    script:
    """
    java -Xmx${memory_limit}G -jar /opt/fdrbench/fdrbench.jar -db ${fasta} -o ${fasta.baseName}-entrapment.fasta -fold ${fold} -level protein -entrapment_label ENTRAPMENT_ -entrapment_pos 0 -uniprot -check
    # 'Reheader' to add entrapment index to database and accession part of the header
    # and remove empty entrapment sequences (which can appear if the original sequence has many Xs)
    sed -r -i -e "s;^>ENTRAPMENT_(.+)\\|(.+)\\|(.+)_([0-9]+)\$;>ENTRAPMENT_\\4_\\1|ENTRAPMENT_\\4_\\2|\\3_\\4;g" -e '\$!N;/>.*\\n\$/d;P;D'  ${fasta.baseName}-entrapment.fasta
    """
}
