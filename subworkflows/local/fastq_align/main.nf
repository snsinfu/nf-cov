/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads with the selected aligner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_BWA     } from '../align_bwa'
include { ALIGN_BWAMEM2 } from '../align_bwamem2'
include { ALIGN_BWAMEM3 } from '../align_bwamem3'

workflow FASTQ_ALIGN {
    take:
    reads    // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    index    // channel: [ val(meta), path(index) ]
    fasta    // channel: path(fasta)
    aligner  // string: 'bwa' | 'bwa-mem2' | 'bwa-mem3'

    main:
    ch_bam = channel.empty()

    // Process/workflow references are not callables, so dispatch with explicit if/else.
    if (aligner == 'bwa') {
        ALIGN_BWA(reads, index, fasta)
        ch_bam = ALIGN_BWA.out.bam
    }
    else if (aligner == 'bwa-mem3') {
        ALIGN_BWAMEM3(reads, index, fasta)
        ch_bam = ALIGN_BWAMEM3.out.bam
    }
    else {
        ALIGN_BWAMEM2(reads, index, fasta)
        ch_bam = ALIGN_BWAMEM2.out.bam
    }

    emit:
    bam = ch_bam   // channel: [ val(meta), path(*.bam) ]
}
