nextflow.enable.dsl=2

percolator_image = 'quay.io/medbioinf/percolator:3.6.5'

// number of threads used by percolator
params.percolator_threads = 4

/**
 * Executes percolator for the given PIN files
 *
 * @return percolated PSMs (TSV pout file)
 */
workflow psm_percolator {
    take:
        pin_files

    main:
        pout_files = run_percolator(pin_files)

    emit:
        pout_files
}


process run_percolator {
    cpus  { params.percolator_threads }
    memory '8GB'

    container { percolator_image }

    input:
    path pin_file

    output:
    path "${pin_file.baseName}.pout"

    """
    percolator --num-threads ${params.percolator_threads} --only-psms --results-psms ${pin_file.baseName}.pout ${pin_file} 
    """
}
