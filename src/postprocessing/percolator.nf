nextflow.enable.dsl=2

params.percolator_image = 'ghcr.io/percolator/percolator:branch-3-08'

// number of threads used by percolator
params.percolator_threads = 4
params.percolator_mem = "4 GB"

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
    memory { params.percolator_mem }
    container { params.percolator_image }

    input:
    path pin_file

    output:
    path "${pin_file.baseName}.pout"

    script:
    """
    percolator --num-threads ${params.percolator_threads} --only-psms --post-processing-tdc --search-input concatenated --results-psms ${pin_file.baseName}.pout ${pin_file} 
    """
}
