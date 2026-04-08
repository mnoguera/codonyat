process BOWTIE2_BUILD {
    label 'process_medium'

    input:
    path(reference)

    output:
    path("bowtie2_index"), emit: index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build ${reference} bowtie2_index/reference
    """
}

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log_out

    script:
    """
    bowtie2 \\
        ${params.bowtie2_args} \\
        -p ${task.cpus} \\
        -x ${index}/reference \\
        ${reads.size() > 1 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"} \\
        2> ${meta.id}_bowtie2.log \\
        | samtools view -@ ${task.cpus} -bS - \\
        > ${meta.id}.bam
    """
}
