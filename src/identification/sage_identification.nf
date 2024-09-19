nextflow.enable.dsl=2

sage_image = 'quay.io/medbioinf/sage:v0.14.7'
python_image = 'medbioinf/ident-comparison-python'

// number of threads used by sage
params.sage_threads = 16

include {convert_and_enhance_psm_tsv} from workflow.projectDir + '/src/postprocessing/convert_and_enhance_psm_tsv.nf'
include {target_decoy_approach} from workflow.projectDir + '/src/postprocessing/default_target_decoy_approach.nf'
include {psm_percolator} from workflow.projectDir + '/src/postprocessing/percolator.nf'
include {ms2rescore_workflow} from workflow.projectDir + '/src/postprocessing/ms2rescore.nf'

/**
 * Executes the identification using Sage
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow sage_identification {
    take:
        sage_config_file
        fasta
        mzmls

    main:
        // all mzMLs are processed at once in sage for now -> much faster
        // this will break as soon as in one SDRF different params are given
        sage_results = identification_with_sage(sage_config_file, fasta, mzmls.collect())
        separated_results = separate_sage_results(sage_results.sage_pin, sage_results.sage_tsv)

        // transpose to tuples containing [pin, tsv] files for each mzML 
        return_files = separated_results.sage_pin.collect()
            .concat(separated_results.sage_tsv.collect())
            .toList()
            .transpose()
        
        psm_tsvs_and_pin = convert_and_enhance_psm_tsv(separated_results.sage_tsv.flatten(), 'sage_tsv', 'sage')
        psm_tsvs = psm_tsvs_and_pin[0]
        pin_files = psm_tsvs_and_pin[1]

        tda_results = target_decoy_approach(psm_tsvs, 'sage')

        pout_files = psm_percolator(pin_files)

        psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.sage')) ] }
        ms2rescore_results = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), 'sage')
        mokapot_results = ms2rescore_results[0]
        mokapot_features = ms2rescore_results[1]
        
    emit:
        return_files
        psm_tsvs
        tda_results
        pout_files
        mokapot_results
        mokapot_features
}


process identification_with_sage {
    cpus { params.sage_threads }
    container { sage_image }

    input:
    path sage_config_file
    path fasta
    path mzmls

    output:
    path "results.sage.pin", emit: sage_pin
    path "results.sage.tsv", emit: sage_tsv

    """
    RAYON_NUM_THREADS=${params.sage_threads} sage ${sage_config_file} -f ${fasta} --write-pin --output_directory ./ ${mzmls}
    """
}


process separate_sage_results {
    container { python_image }

    input:
    path sage_pin
    path sage_tsv

    output:
    path "*.sage.pin", emit: sage_pin
    path "*.sage.tsv", emit: sage_tsv

    """
    # process the pin file and create one file for each input file
    for filename in \$(awk 'NR>1{a[\$6]++} END{for(b in a) print b}' ${sage_pin});
    do
        head -n1 ${sage_pin} > \${filename}.sage.pin
        awk -v f1="\${filename}" '\$6==f1' ${sage_pin} >> \${filename}.sage.pin
    done

    # process the tsv file and create one file for each input file
    for filename in \$(awk 'NR>1{a[\$5]++} END{for(b in a) print b}' ${sage_tsv});
    do
        head -n1 ${sage_tsv} > \${filename}.sage.tsv
        awk -v f1="\${filename}" '\$5==f1' ${sage_tsv} >> \${filename}.sage.tsv
    done
    """
}
