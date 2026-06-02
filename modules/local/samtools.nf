process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/env/samtools.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_sorted.bam"), emit: bam

    script:
    """
    # Validate input
    if [ ! -f ${bam} ]; then
        echo "ERROR: Input BAM not found: ${bam}" >&2
        exit 1
    fi

    # Sort BAM
    samtools sort -@ ${task.cpus} -o ${meta.id}_sorted.bam ${bam} 2> sort.log || {
        echo "ERROR: samtools sort failed for ${meta.id}" >&2
        cat sort.log >&2
        exit 1
    }

    # Validate output
    if [ ! -f ${meta.id}_sorted.bam ]; then
        echo "ERROR: Sorting produced no output file" >&2
        exit 1
    fi
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/env/samtools.yml"

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
    conda "${projectDir}/env/samtools.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    """
    # Validate input
    if [ ! -f ${bam} ]; then
        echo "ERROR: Input BAM not found: ${bam}" >&2
        exit 1
    fi

    # Convert to SAM
    samtools view -h -o ${meta.id}.sam ${bam} 2> view.log || {
        echo "ERROR: BAM to SAM conversion failed for ${meta.id}" >&2
        cat view.log >&2
        exit 1
    }

    # Validate output is non-empty
    if [ ! -s ${meta.id}.sam ]; then
        echo "ERROR: Conversion produced empty SAM file" >&2
        exit 1
    fi
    """
}
