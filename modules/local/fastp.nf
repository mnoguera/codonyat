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
    fastp \\
        -i ${reads[0]} \\
        ${reads.size() > 1 ? "-I ${reads[1]}" : ""} \\
        -o ${meta.id}_R1_trimmed.fastq.gz \\
        ${reads.size() > 1 ? "-O ${meta.id}_R2_trimmed.fastq.gz" : ""} \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        --thread ${task.cpus} \\
        ${params.fastp_args} \\
        2> ${meta.id}_fastp.log
    """
}
