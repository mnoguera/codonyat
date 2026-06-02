class Utils {

    static void validateSamplesheet(String samplesheetPath) {
        def samplesheet = new File(samplesheetPath)

        if (!samplesheet.exists()) {
            throw new Exception("Samplesheet not found: ${samplesheetPath}")
        }

        def lines = samplesheet.readLines()
        if (lines.size() < 2) {
            throw new Exception("Samplesheet is empty or has no data rows")
        }

        // Note: Simple split(',') doesn't handle quoted CSV values with commas inside them
        // For complex CSV parsing, consider using a dedicated CSV library
        def header = lines[0].split(',')
        def requiredCols = ['sample_id', 'fastq_1', 'fastq_2']

        requiredCols.each { col ->
            if (!header.contains(col)) {
                throw new Exception("Missing required column in samplesheet: ${col}")
            }
        }

        def sampleIds = [] as Set
        def sheetDir = samplesheet.getParentFile()

        lines.drop(1).eachWithIndex { line, idx ->
            def fields = line.split(',')

            // Validate field count before accessing array elements
            if (fields.size() < header.size()) {
                throw new Exception("Row ${idx + 2}: insufficient fields (expected ${header.size()}, got ${fields.size()})")
            }

            def sampleId = fields[header.indexOf('sample_id')]
            def fastq1 = fields[header.indexOf('fastq_1')]
            def fastq2 = fields[header.indexOf('fastq_2')]

            // Check for duplicate sample IDs
            if (sampleIds.contains(sampleId)) {
                throw new Exception("Row ${idx + 2}: Duplicate sample_id in samplesheet: ${sampleId}")
            }
            sampleIds.add(sampleId)

            // Validate fastq_1 exists
            def r1File = resolvePath(fastq1, sheetDir)
            if (!r1File.exists()) {
                throw new Exception("Row ${idx + 2}: FASTQ file not found: ${fastq1} (sample: ${sampleId})")
            }
            if (!isValidFastqExtension(fastq1)) {
                throw new Exception("Row ${idx + 2}: Invalid FASTQ extension: ${fastq1} (must be .fastq, .fq, .fastq.gz, .fq.gz)")
            }

            // Validate fastq_2 if paired-end
            if (fastq2 && fastq2.trim()) {
                def r2File = resolvePath(fastq2, sheetDir)
                if (!r2File.exists()) {
                    throw new Exception("Row ${idx + 2}: FASTQ file not found: ${fastq2} (sample: ${sampleId})")
                }
                if (!isValidFastqExtension(fastq2)) {
                    throw new Exception("Row ${idx + 2}: Invalid FASTQ extension: ${fastq2}")
                }
            }
        }
    }

    static void validateInputFile(String filePath, String fileType) {
        def file = new File(filePath)
        if (!file.exists()) {
            throw new Exception("${fileType} file not found: ${filePath}")
        }
        if (!file.canRead()) {
            throw new Exception("${fileType} file is not readable: ${filePath}")
        }
    }

    static File resolvePath(String path, File baseDir) {
        def file = new File(path)
        return file.isAbsolute() ? file : new File(baseDir, path)
    }

    static boolean isValidFastqExtension(String filename) {
        return filename =~ /\.(fastq|fq)(\.gz)?$/
    }
}
