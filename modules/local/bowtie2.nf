process BOWTIE2_BUILD {
    label 'process_medium'
    conda "${projectDir}/env/bowtie2.yml"

    input:
    path(reference)

    output:
    path("bowtie2_index"), emit: index

    script:
    """
    # Validate input
    if [ ! -f ${reference} ]; then
        echo "ERROR: Reference file not found: ${reference}" >&2
        exit 1
    fi

    # Build index
    mkdir -p bowtie2_index
    bowtie2-build ${reference} bowtie2_index/reference 2> build.log || {
        echo "ERROR: bowtie2-build failed. Invalid FASTA format?" >&2
        cat build.log >&2
        exit 1
    }

    # Validate output
    if [ ! -f bowtie2_index/reference.1.bt2 ]; then
        echo "ERROR: Index build produced no output files" >&2
        exit 1
    fi
    """
}

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/env/bowtie2.yml"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log_out

    script:
    """
    # Validate index exists
    if [ ! -f ${index}/reference.1.bt2 ]; then
        echo "ERROR: Bowtie2 index not found in ${index}" >&2
        exit 1
    fi

    # Run alignment
    bowtie2 \\
        ${params.bowtie2_args} \\
        -p ${task.cpus} \\
        -x ${index}/reference \\
        ${reads.size() > 1 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"} \\
        2> ${meta.id}_bowtie2.log \\
        | samtools view -@ ${task.cpus} -bS - \\
        > ${meta.id}.bam || {
        echo "ERROR: Bowtie2 alignment failed for ${meta.id}" >&2
        cat ${meta.id}_bowtie2.log >&2
        exit 1
    }

    # Check for alignments (warning only, not fatal)
    NUM_ALIGNMENTS=\$(samtools view -c ${meta.id}.bam)
    if [ "\$NUM_ALIGNMENTS" -eq 0 ]; then
        echo "WARNING: Zero alignments for ${meta.id}. Wrong reference or corrupted reads?" >&2
    fi
    """
}
