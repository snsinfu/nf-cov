/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome-wide coverage bigWigs with deepTools bamCoverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEEPTOOLS_BAMCOVERAGE } from '../../../modules/nf-core/deeptools/bamcoverage'

workflow BAM_COVERAGE {
    take:
    ch_bam    // channel: [ val(meta), path(*.bam), path(*.bai) ]
    ch_ref    // channel: [ val(egs), path(fasta), path(fai) ] (single element)

    main:
    //
    // Broadcast the effective genome size / reference across the BAM queue.
    //
    ch_bam_egs = ch_bam
        .combine(ch_ref)
        .map { meta, bam, bai, egs, fasta, fai -> [ meta + [ egs: egs, fasta: fasta, fai: fai ], bam, bai ] }

    DEEPTOOLS_BAMCOVERAGE(
        ch_bam_egs.map { meta, bam, bai -> [ meta, bam, bai ] },
        ch_bam_egs.map { meta, _bam, _bai -> meta.fasta }.first(),
        ch_bam_egs.map { meta, _bam, _bai -> meta.fai }.first(),
        [ [], [] ]
    )

    emit:
    bigwig = DEEPTOOLS_BAMCOVERAGE.out.bigwig   // channel: [ val(meta), path(*.bigWig) ]
}
