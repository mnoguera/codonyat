process BOWTIE2_BUILD {
    label 'process_medium'
    container 'biocontainers/bowtie2:2.5.4--he20e202_1'

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
    container 'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f7f8571f6e1e50f8f-0'

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
