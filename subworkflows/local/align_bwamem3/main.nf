/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads with bwa-mem3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM3_MEM } from '../../../modules/nf-core/bwamem3/mem'

workflow ALIGN_BWAMEM3 {
    take:
    reads   // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    index   // channel: [ val(meta), path(index) ]
    fasta   // channel: path(fasta)

    main:
    BWAMEM3_MEM(reads, index, fasta.map { f -> [ [:], f ] }, true)

    emit:
    // nf-core bwamem3/mem emits the alignment under `aligned`, not `bam`
    bam = BWAMEM3_MEM.out.aligned   // channel: [ val(meta), path(*.bam) ]
}
