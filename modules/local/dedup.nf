process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_dedup.bam"),         emit: bam
    tuple val(meta), path("*_dedup_metrics.txt"), emit: metrics

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${meta.id}_dedup.bam \\
        METRICS_FILE=${meta.id}_dedup_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}
