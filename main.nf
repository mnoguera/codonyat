#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Validate required params
if (!params.input)     { error "Please provide --input samplesheet.csv" }
if (!params.reference) { error "Please provide --reference reference.fasta" }
if (!params.amplicons) { error "Please provide --amplicons amplicons.tsv" }

// Include modules
include { FASTQC as FASTQC_RAW     } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/local/fastqc'
include { FASTP                     } from './modules/local/fastp'
include { BOWTIE2_BUILD             } from './modules/local/bowtie2'
include { BOWTIE2_ALIGN             } from './modules/local/bowtie2'
include { SAMTOOLS_SORT             } from './modules/local/samtools'
include { BAM_TO_SAM                } from './modules/local/samtools'
include { PICARD_MARKDUPLICATES     } from './modules/local/dedup'
include { CODONYAT                  } from './modules/local/codonyat'
include { MULTIQC                   } from './modules/local/multiqc'

workflow {
    // Parse samplesheet — resolve relative FASTQ paths against the samplesheet's directory
    def sheet = file(params.input)
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id]
            def r1   = file(row.fastq_1.startsWith('/') ? row.fastq_1 : "${sheet.parent}/${row.fastq_1}", checkIfExists: true)
            def reads = row.fastq_2 ? [r1, file(row.fastq_2.startsWith('/') ? row.fastq_2 : "${sheet.parent}/${row.fastq_2}", checkIfExists: true)] : [r1]
            [meta, reads]
        }

    ch_reference = file(params.reference, checkIfExists: true)
    ch_amplicons = file(params.amplicons, checkIfExists: true)

    // Collect QC outputs for MultiQC
    ch_multiqc = Channel.empty()

    // 1. Raw read QC
    if (!params.skip_fastqc) {
        FASTQC_RAW(ch_input)
        ch_multiqc = ch_multiqc.mix(FASTQC_RAW.out.zip.map{ meta, zip -> zip })
    }

    // 2. Trim and filter
    FASTP(ch_input)
    ch_multiqc = ch_multiqc.mix(FASTP.out.json.map{ meta, json -> json })

    // 3. Trimmed read QC
    if (!params.skip_fastqc) {
        FASTQC_TRIMMED(FASTP.out.reads)
        ch_multiqc = ch_multiqc.mix(FASTQC_TRIMMED.out.zip.map{ meta, zip -> zip })
    }

    // 4. Build bowtie2 index
    BOWTIE2_BUILD(ch_reference)

    // 5. Align (--no-unal filters non-reference reads = decontamination)
    BOWTIE2_ALIGN(FASTP.out.reads, BOWTIE2_BUILD.out.index)
    ch_multiqc = ch_multiqc.mix(BOWTIE2_ALIGN.out.log_out.map{ meta, log -> log })

    // 6. Sort
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)

    // 7. Optional deduplication
    if (params.skip_dedup) {
        ch_final_bam = SAMTOOLS_SORT.out.bam
    } else {
        PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)
        ch_final_bam = PICARD_MARKDUPLICATES.out.bam
        ch_multiqc = ch_multiqc.mix(PICARD_MARKDUPLICATES.out.metrics.map{ meta, m -> m })
    }

    // 8. Convert BAM to SAM for codonyat
    BAM_TO_SAM(ch_final_bam)

    // 9. Variant calling
    CODONYAT(BAM_TO_SAM.out.sam, ch_reference, ch_amplicons)

    // 10. MultiQC
    MULTIQC(ch_multiqc.collect().ifEmpty([]))
}
