process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),             emit: json
    tuple val(meta), path("*.html"),             emit: html
    tuple val(meta), path("*.log"),              emit: log_out

    script:
    """
    # Validate input exists
    if [ ! -f ${reads[0]} ]; then
        echo "ERROR: Input FASTQ not found: ${reads[0]}" >&2
        exit 1
    fi

    # Run fastp
    fastp \\
        -i ${reads[0]} \\
        ${reads.size() > 1 ? "-I ${reads[1]}" : ""} \\
        -o ${meta.id}_R1_trimmed.fastq.gz \\
        ${reads.size() > 1 ? "-O ${meta.id}_R2_trimmed.fastq.gz" : ""} \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        --thread ${task.cpus} \\
        ${params.fastp_args} \\
        2> ${meta.id}_fastp.log || {
        echo "ERROR: fastp failed for ${meta.id}" >&2
        cat ${meta.id}_fastp.log >&2
        exit 1
    }

    # Validate output exists and is non-empty
    if [ ! -s ${meta.id}_R1_trimmed.fastq.gz ]; then
        echo "ERROR: Trimming produced no output or empty file" >&2
        exit 1
    fi
    """
}
