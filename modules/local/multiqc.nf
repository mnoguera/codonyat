process MULTIQC {
    label 'process_low'

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
