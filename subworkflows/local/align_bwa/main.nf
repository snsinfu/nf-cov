/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads with BWA-MEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_MEM } from '../../../modules/nf-core/bwa/mem'

workflow ALIGN_BWA {
    take:
    reads   // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    index   // channel: [ val(meta), path(index) ]
    fasta   // channel: path(fasta)

    main:
    BWA_MEM(reads, index, fasta.map { f -> [ [:], f ] }, true)

    emit:
    bam = BWA_MEM.out.bam   // channel: [ val(meta), path(*.bam) ]
}
