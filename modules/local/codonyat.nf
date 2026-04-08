process CODONYAT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sam)
    path(reference)
    path(amplicons)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.xml"), emit: xml

    script:
    """
    codonyat ${sam} ${reference} ${amplicons} \\
        --protein ${params.protein} \\
        --ratio-upper ${params.ratio_upper} \\
        --ratio-lower ${params.ratio_lower} \\
        --entropy-threshold ${params.entropy_threshold}
    """
}
