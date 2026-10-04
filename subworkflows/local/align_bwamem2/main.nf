/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads with bwa-mem2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM } from '../../../modules/nf-core/bwamem2/mem'

workflow ALIGN_BWAMEM2 {
    take:
    reads   // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    index   // channel: [ val(meta), path(index) ]
    fasta   // channel: path(fasta)

    main:
    BWAMEM2_MEM(reads, index, fasta.map { f -> [ [:], f ] }, true)

    emit:
    bam = BWAMEM2_MEM.out.bam   // channel: [ val(meta), path(*.bam) ]
}
