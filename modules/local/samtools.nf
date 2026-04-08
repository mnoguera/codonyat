process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    container 'biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}_sorted.bam ${bam}
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'
    container 'biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai

    script:
    """
    samtools index ${bam}
    """
}

process BAM_TO_SAM {
    tag "${meta.id}"
    label 'process_low'
    container 'biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    """
    samtools view -h -o ${meta.id}.sam ${bam}
    """
}
