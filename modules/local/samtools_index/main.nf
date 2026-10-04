/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Index a BAM/SAM/CRAM file

    Local fork of nf-core/modules/samtools/index (@ efec54255f9baad3ea032173d75031929883bed8)
    that adds a `val index_format` input. Upstream has no such input and selects the format
    only through `ext.args`; the explicit input keeps the two index producers in nf-cov
    (this process and samtools/view) on one symmetric interface.
    Remove this fork if upstream gains an `index_format` input.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(input)
    val index_format // 'bai' | 'csi' | 'crai'

    output:
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def format = index_format ?: 'bai'
    if (!format.matches('bai|csi|crai')) {
        error("Index format not one of bai, csi, crai.")
    }
    def format_arg = format == 'csi' ? '-c' : format == 'bai' ? '-b' : ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${format_arg} \\
        ${args} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def format = index_format ?: 'bai'
    if (!format.matches('bai|csi|crai')) {
        error("Index format not one of bai, csi, crai.")
    }
    def extension = file(input).getExtension() == 'cram' ? "crai" : format
    """
    touch ${input}.${extension}
    """
}
