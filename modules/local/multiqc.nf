process MULTIQC {
    label 'process_low'
    container 'multiqc/multiqc:1.22.2'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    multiqc . --filename multiqc_report
    """
}
