process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
