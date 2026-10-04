/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'

workflow NFCOV {
    take:
    fasta         // path: genome fasta
    index         // path: prebuilt aligner index (optional)
    aligner       // string: 'bwa' | 'bwa-mem2' | 'bwa-mem3'
    catalog_gsize // val: map read-length -> gsize, or null
    explicit_egs  // val: int, or null
    read_length   // val: int, or null

    main:
    //
    // Reference preparation
    //
    PREPARE_GENOME(fasta, index, aligner, catalog_gsize, explicit_egs, read_length)

    emit:
    fasta = PREPARE_GENOME.out.fasta   // channel: path(genome.fa)
    fai   = PREPARE_GENOME.out.fai     // channel: [ val(meta), path(genome.fa.fai) ]
    index = PREPARE_GENOME.out.index   // channel: [ val(meta), path(index directory) ]
    egs   = PREPARE_GENOME.out.egs     // channel: val(effective genome size)
}
