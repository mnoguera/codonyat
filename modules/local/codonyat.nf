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
    # Validate inputs
    if [ ! -s ${sam} ]; then
        echo "ERROR: SAM file is empty or missing: ${sam}" >&2
        exit 1
    fi

    if [ ! -f ${reference} ]; then
        echo "ERROR: Reference file not found: ${reference}" >&2
        exit 1
    fi

    if [ ! -f ${amplicons} ]; then
        echo "ERROR: Amplicons file not found: ${amplicons}" >&2
        exit 1
    fi

    # Run codonyat with error capture
    codonyat ${sam} ${reference} ${amplicons} \\
        --protein ${params.protein} \\
        --ratio-upper ${params.ratio_upper} \\
        --ratio-lower ${params.ratio_lower} \\
        --entropy-threshold ${params.entropy_threshold} 2>&1 | tee codonyat.log

    EXIT_CODE=\${PIPESTATUS[0]}
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "ERROR: codonyat failed with exit code \$EXIT_CODE" >&2
        echo "Common causes: invalid protein name, mismatched amplicon coordinates" >&2
        cat codonyat.log >&2
        exit 1
    fi

    # Validate outputs were created
    TSV_COUNT=\$(ls *.tsv 2>/dev/null | wc -l)
    XML_COUNT=\$(ls *.xml 2>/dev/null | wc -l)
    if [ "\$TSV_COUNT" -eq 0 ] || [ "\$XML_COUNT" -eq 0 ]; then
        echo "ERROR: codonyat did not produce expected output files" >&2
        exit 1
    fi
    """
}
