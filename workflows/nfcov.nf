/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCOV {

    main:
    ch_versions = channel.empty()

    // Pipeline steps are wired incrementally (see .opencode/plans/).

    emit:
    versions = ch_versions   // channel: [ path(versions.yml) ]
}
