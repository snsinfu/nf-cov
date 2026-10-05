#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ to genome-wide coverage bigWig.
----------------------------------------------------------------------------------------
*/

include { NFCOV } from './workflows/nfcov'

workflow {
    //
    // Resolve reference paths from the genome catalog only when not explicitly given, so
    // CLI/config --fasta/--bwa_index/--bwamem2_index/--bwamem3_index always win.
    //
    def fasta = params.fasta ?: getGenomeAttribute('fasta')
    if (!fasta) {
        error("Genome fasta file not specified: use --genome <key>, --fasta <file.fa>, or a custom config.")
    }
    if (!(params.aligner in ['bwa', 'bwa-mem2', 'bwa-mem3'])) {
        error("Invalid --aligner '${params.aligner}'. Use 'bwa', 'bwa-mem2' or 'bwa-mem3'.")
    }

    //
    // Validate the paired-end fragment-size filter bounds (null disables that direction).
    //
    def min_frag = params.min_fragment_size != null ? params.min_fragment_size as Integer : null
    def max_frag = params.max_fragment_size != null ? params.max_fragment_size as Integer : null
    if (min_frag != null && min_frag < 0) {
        error("--min_fragment_size must be >= 0 (got ${min_frag}).")
    }
    if (max_frag != null && max_frag < 0) {
        error("--max_fragment_size must be >= 0 (got ${max_frag}).")
    }
    if (min_frag != null && max_frag != null && min_frag > max_frag) {
        error("--min_fragment_size (${min_frag}) must be <= --max_fragment_size (${max_frag}).")
    }

    def index = resolveAlignerIndex(params.aligner)

    //
    // Warn when another aligner's index parameter is supplied but ignored.
    //
    def index_params = [
        'bwa'      : [flag: '--bwa_index',     value: params.bwa_index],
        'bwa-mem2' : [flag: '--bwamem2_index', value: params.bwamem2_index],
        'bwa-mem3' : [flag: '--bwamem3_index', value: params.bwamem3_index]
    ]
    def own_index = index_params[params.aligner]
    if (own_index && !own_index.value) {
        index_params.each { name, spec ->
            if (name != params.aligner && spec.value) {
                log.warn("${spec.flag} is ignored with --aligner ${params.aligner} (bwa, bwa-mem2 and bwa-mem3 indexes are not interchangeable). Use ${own_index.flag} for a precomputed ${params.aligner} index.")
            }
        }
    }

    def catalog_gsize = getGenomeAttribute('macs_gsize')

    NFCOV(
        fasta,
        index,
        params.aligner,
        catalog_gsize,
        params.effective_genome_size,
        params.read_length
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get an attribute from the genome catalog, e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Resolve the precomputed index for an aligner. bwa, bwa-mem2 and bwa-mem3 indexes are not
// interchangeable, so each aligner reads only its own parameter (falling back to its own
// catalog key).
//
def resolveAlignerIndex(aligner) {
    if (aligner == 'bwa') {
        return params.bwa_index     ?: getGenomeAttribute('bwa')
    }
    else if (aligner == 'bwa-mem2') {
        return params.bwamem2_index ?: getGenomeAttribute('bwamem2')
    }
    else if (aligner == 'bwa-mem3') {
        return params.bwamem3_index ?: getGenomeAttribute('bwamem3')
    }
    return null
}
