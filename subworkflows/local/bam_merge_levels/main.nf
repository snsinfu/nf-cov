/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Replicate merges:
      .mLb  merge technical replicates within a biological replicate
      .mRp  merge biological replicates of a sample (only when >1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_MLB } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_MRP } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MLB } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MRP } from '../../../modules/nf-core/samtools/index'

workflow BAM_MERGE_LEVELS {
    take:
    ch_bam        // channel: [ val(meta library), path(clean *.bam) ]
    ch_ref        // channel: [ path(fasta), path(fai) ] (single element)
    _save_library // val: boolean (publishing handled via config)

    main:
    ch_merge_ref = ch_ref.map { f, fi -> [ [:], f, fi, [] ] }

    //
    // .mLb: merge technical replicates of a (sample, biological replicate)
    //
    ch_mlb_in = ch_bam
        .map { meta, bam -> [ "${meta.sample}_REP${meta.replicate}".toString(), meta, bam ] }
        .groupTuple(by: [0])
        .map { mlb_id, metas, bams ->
            def m = [
                id         : mlb_id,
                sample     : metas[0].sample,
                replicate  : metas[0].replicate,
                single_end : metas[0].single_end,
                level      : 'mLb'
            ]
            [ m, bams.flatten(), [] ]
        }

    SAMTOOLS_MERGE_MLB(ch_mlb_in, ch_merge_ref)
    ch_mlb = SAMTOOLS_MERGE_MLB.out.bam
        .map { meta, bam ->
            def m = meta + [ id: "${meta.id}.mLb".toString(), library: meta.id ]
            [ m, bam ]
        }
    SAMTOOLS_INDEX_MLB(ch_mlb)
    ch_mlb_bam_bai = ch_mlb.join(SAMTOOLS_INDEX_MLB.out.index, by: [0])

    //
    // .mRp: merge biological replicates of a sample, only when there is more than one
    //
    ch_mrp = channel.empty()
    ch_mrp_in = ch_mlb
        .map { meta, bam -> [ meta.sample.toString(), meta, bam ] }
        .groupTuple(by: [0])
        .filter { _sample, _metas, bams -> bams.flatten().size() > 1 }
        .map { sample, metas, bams ->
            def m = [
                id         : "${sample}.mRp".toString(),
                sample     : sample,
                single_end : metas[0].single_end,
                level      : 'mRp'
            ]
            [ m, bams.flatten(), [] ]
        }

    SAMTOOLS_MERGE_MRP(ch_mrp_in, ch_merge_ref)
    ch_mrp = SAMTOOLS_MERGE_MRP.out.bam
    SAMTOOLS_INDEX_MRP(ch_mrp)
    ch_mrp_bam_bai = ch_mrp.join(SAMTOOLS_INDEX_MRP.out.index, by: [0])

    emit:
    mLb = ch_mlb_bam_bai   // channel: [ val(meta), path(*.mLb.bam), path(*.bai) ]
    mRp = ch_mrp_bam_bai   // channel: [ val(meta), path(*.mRp.bam), path(*.bai) ]
}
