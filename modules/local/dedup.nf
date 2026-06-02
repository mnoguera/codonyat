process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/env/picard.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_dedup.bam"),         emit: bam
    tuple val(meta), path("*_dedup_metrics.txt"), emit: metrics

    script:
    """
    # Validate input
    if [ ! -f ${bam} ]; then
        echo "ERROR: Input BAM not found: ${bam}" >&2
        exit 1
    fi

    # Run Picard MarkDuplicates
    java -jar /usr/picard/picard.jar MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${meta.id}_dedup.bam \\
        METRICS_FILE=${meta.id}_dedup_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        VALIDATION_STRINGENCY=LENIENT 2> picard.log || {
        echo "ERROR: Picard MarkDuplicates failed for ${meta.id}" >&2
        cat picard.log >&2
        exit 1
    }

    # Validate outputs
    if [ ! -f ${meta.id}_dedup.bam ]; then
        echo "ERROR: Deduplication produced no BAM output" >&2
        exit 1
    fi

    if [ ! -f ${meta.id}_dedup_metrics.txt ]; then
        echo "ERROR: Deduplication produced no metrics file" >&2
        exit 1
    fi
    """
}
