process MULTIQC {
    label 'process_low'
    conda "${projectDir}/env/multiqc.yml"

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    """
    multiqc . --filename multiqc_report
    """
}
