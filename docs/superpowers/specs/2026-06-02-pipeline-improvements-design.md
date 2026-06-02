# Pipeline Reliability and Maintainability Improvements

**Date:** 2026-06-02  
**Author:** Claude Code  
**Status:** Approved

## Overview

This design addresses five high-impact improvements to the codonyat Nextflow pipeline: retry strategy, input validation, error handling, conda profile support, and Dockerfile improvements. The goal is to enhance pipeline reliability in production while making local development easier.

## Motivation

The current pipeline is functional but has gaps that cause production pain:

1. **Transient failures waste compute time** - Network hiccups or I/O errors kill entire runs with no retry
2. **Late error detection** - Invalid inputs aren't caught until hours into execution
3. **Cryptic failure messages** - When processes fail, logs don't explain why
4. **Docker-only development** - No fast local iteration without containers
5. **Version mismatch risk** - Dockerfile installs PyPI version, not local code

These improvements target the 80/20 wins: maximum reliability and developer experience improvement without over-engineering.

## Design

### 1. Retry Strategy

**Goal:** Automatically retry processes that fail due to transient issues (network, I/O, resource contention).

**Implementation:**

Add retry configuration to `conf/base.config` using Nextflow's `errorStrategy` and `maxRetries` directives:

```groovy
process {
    // High-compute processes: 2 retries (alignment, variant calling)
    withLabel: 'process_high' {
        cpus   = 8
        memory = 16.GB
        time   = 8.h
        errorStrategy = { task.attempt <= 2 ? 'retry' : 'finish' }
        maxRetries = 2
    }
    
    // Medium processes: 3 retries (QC, trimming, format conversion)
    withLabel: 'process_medium' {
        cpus   = 4
        memory = 8.GB
        time   = 4.h
        errorStrategy = { task.attempt <= 3 ? 'retry' : 'finish' }
        maxRetries = 3
    }
    
    // Light processes: 3 retries (metadata operations)
    withLabel: 'process_low' {
        cpus   = 2
        memory = 4.GB
        time   = 1.h
        errorStrategy = { task.attempt <= 3 ? 'retry' : 'finish' }
        maxRetries = 3
    }
}
```

**Rationale:**
- Heavy compute tasks (alignment) get fewer retries because failures are usually deterministic
- I/O-heavy tasks (QC, trimming) get more retries because network/disk issues are often transient
- `errorStrategy` uses closure to check attempt count, falling back to `'finish'` (pipeline continues with other samples)

**Files affected:**
- `conf/base.config` - add retry directives to existing process labels

### 2. Input Validation

**Goal:** Validate samplesheet format and file existence before pipeline execution, catching errors in seconds instead of hours.

**Implementation:**

Create `lib/Utils.groovy` with validation functions:

```groovy
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
            def sampleId = fields[header.indexOf('sample_id')]
            def fastq1 = fields[header.indexOf('fastq_1')]
            def fastq2 = fields[header.indexOf('fastq_2')]
            
            // Check for duplicate sample IDs
            if (sampleIds.contains(sampleId)) {
                throw new Exception("Duplicate sample_id in samplesheet: ${sampleId}")
            }
            sampleIds.add(sampleId)
            
            // Validate fastq_1 exists
            def r1File = resolvePath(fastq1, sheetDir)
            if (!r1File.exists()) {
                throw new Exception("FASTQ file not found: ${fastq1} (sample: ${sampleId})")
            }
            if (!isValidFastqExtension(fastq1)) {
                throw new Exception("Invalid FASTQ extension: ${fastq1} (must be .fastq, .fq, .fastq.gz, .fq.gz)")
            }
            
            // Validate fastq_2 if paired-end
            if (fastq2 && fastq2.trim()) {
                def r2File = resolvePath(fastq2, sheetDir)
                if (!r2File.exists()) {
                    throw new Exception("FASTQ file not found: ${fastq2} (sample: ${sampleId})")
                }
                if (!isValidFastqExtension(fastq2)) {
                    throw new Exception("Invalid FASTQ extension: ${fastq2}")
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
```

Call validation in `main.nf` before channel creation:

```groovy
workflow {
    // Validate inputs before starting
    Utils.validateSamplesheet(params.input)
    Utils.validateInputFile(params.reference, "Reference FASTA")
    Utils.validateInputFile(params.amplicons, "Amplicons TSV")
    
    // Parse samplesheet...
    // (existing code)
}
```

**Validation checks:**
- Samplesheet exists and has data
- Required columns present
- No duplicate sample IDs
- All FASTQ files exist and have valid extensions
- Paired-end samples have both R1 and R2
- Reference and amplicons files exist and are readable

**Files affected:**
- `lib/Utils.groovy` (new file)
- `main.nf` - add validation calls at workflow start

### 3. Error Handling

**Goal:** Add explicit error messages in process scripts so failures provide actionable debugging information.

**Implementation:**

Add defensive checks and clear error messages to each process module. Pattern:

```bash
# 1. Validate inputs exist and are non-empty
if [ ! -s input_file ]; then
    echo "ERROR: Input file is empty or missing: input_file" >&2
    exit 1
fi

# 2. Run tool with error check
tool --args input_file > output_file 2> tool.log || {
    echo "ERROR: Tool failed. Check tool.log for details." >&2
    cat tool.log >&2
    exit 1
}

# 3. Validate output was created
if [ ! -f output_file ]; then
    echo "ERROR: Tool did not produce expected output: output_file" >&2
    exit 1
fi
```

**Process-specific error handling:**

**BOWTIE2_BUILD:**
```bash
if [ ! -f ${reference} ]; then
    echo "ERROR: Reference file not found: ${reference}" >&2
    exit 1
fi

mkdir -p bowtie2_index
bowtie2-build ${reference} bowtie2_index/reference 2> build.log || {
    echo "ERROR: bowtie2-build failed. Invalid FASTA format?" >&2
    cat build.log >&2
    exit 1
}

if [ ! -f bowtie2_index/reference.1.bt2 ]; then
    echo "ERROR: Index build produced no output files" >&2
    exit 1
fi
```

**BOWTIE2_ALIGN:**
```bash
# Check index exists
if [ ! -f ${index}/reference.1.bt2 ]; then
    echo "ERROR: Bowtie2 index not found in ${index}" >&2
    exit 1
fi

# Run alignment
bowtie2 ${params.bowtie2_args} ... 2> ${meta.id}_bowtie2.log || {
    echo "ERROR: Bowtie2 alignment failed for ${meta.id}" >&2
    cat ${meta.id}_bowtie2.log >&2
    exit 1
}

# Check for alignments
NUM_ALIGNMENTS=\$(samtools view -c ${meta.id}.bam)
if [ "\$NUM_ALIGNMENTS" -eq 0 ]; then
    echo "WARNING: Zero alignments for ${meta.id}. Wrong reference or corrupted reads?" >&2
fi
```

**CODONYAT:**
```bash
# Validate inputs
if [ ! -s ${sam} ]; then
    echo "ERROR: SAM file is empty: ${sam}" >&2
    exit 1
fi

if [ ! -f ${reference} ]; then
    echo "ERROR: Reference file not found: ${reference}" >&2
    exit 1
fi

if [ ! -f ${amplicons} ]; then
    echo "ERROR: Amplicons file not found: ${amplicons}" >&2
    exit 1
fi

# Run codonyat with error capture
codonyat ${sam} ${reference} ${amplicons} \
    --protein ${params.protein} \
    --ratio-upper ${params.ratio_upper} \
    --ratio-lower ${params.ratio_lower} \
    --entropy-threshold ${params.entropy_threshold} 2>&1 | tee codonyat.log

EXIT_CODE=\${PIPESTATUS[0]}
if [ \$EXIT_CODE -ne 0 ]; then
    echo "ERROR: codonyat failed with exit code \$EXIT_CODE" >&2
    echo "Common causes: invalid protein name, mismatched amplicon coordinates" >&2
    cat codonyat.log >&2
    exit 1
fi

# Validate outputs were created
if [ ! -f *.tsv ] || [ ! -f *.xml ]; then
    echo "ERROR: codonyat did not produce expected output files" >&2
    exit 1
fi
```

**FASTP:**
```bash
# Check input exists
if [ ! -f ${reads[0]} ]; then
    echo "ERROR: Input FASTQ not found: ${reads[0]}" >&2
    exit 1
fi

# Run fastp
fastp ... 2> ${meta.id}_fastp.log || {
    echo "ERROR: fastp failed for ${meta.id}" >&2
    cat ${meta.id}_fastp.log >&2
    exit 1
}

# Check trimmed output exists and is non-empty
if [ ! -s ${meta.id}_R1_trimmed.fastq.gz ]; then
    echo "ERROR: Trimming produced no output or empty file" >&2
    exit 1
fi
```

**Files affected:**
- `modules/local/bowtie2.nf`
- `modules/local/codonyat.nf`
- `modules/local/fastp.nf`
- `modules/local/samtools.nf`
- `modules/local/fastqc.nf`
- `modules/local/dedup.nf`

### 4. Conda Profile

**Goal:** Enable local development and execution without Docker using conda/mamba environments.

**Implementation:**

Create `env/` directory with per-tool conda environment YAML files:

**env/fastqc.yml:**
```yaml
name: fastqc
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastqc=0.12.1
```

**env/fastp.yml:**
```yaml
name: fastp
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastp=0.23.4
```

**env/bowtie2.yml:**
```yaml
name: bowtie2
channels:
  - conda-forge
  - bioconda
dependencies:
  - bowtie2=2.5.1
  - samtools=1.19
```

**env/samtools.yml:**
```yaml
name: samtools
channels:
  - conda-forge
  - bioconda
dependencies:
  - samtools=1.19
```

**env/picard.yml:**
```yaml
name: picard
channels:
  - conda-forge
  - bioconda
dependencies:
  - picard=3.1.1
```

**env/codonyat.yml:**
```yaml
name: codonyat
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.10
  - pip
  - biopython>=1.79
  - pip:
    - -e .
```

**env/multiqc.yml:**
```yaml
name: multiqc
channels:
  - conda-forge
  - bioconda
dependencies:
  - multiqc=1.22.2
```

Add conda profile to `nextflow.config`:

```groovy
profiles {
    conda {
        conda.enabled = true
        conda.useMamba = true  // Use mamba for faster solves
    }
    // existing profiles...
}
```

Update each process in modules to specify its conda environment:

```groovy
process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/env/fastqc.yml"
    
    // rest of process...
}

process CODONYAT {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/env/codonyat.yml"
    
    // rest of process...
}
```

**Usage:**
```bash
# Install mamba (faster than conda)
conda install -n base -c conda-forge mamba

# Run with conda
nextflow run . -profile conda \
    --input samplesheet.csv \
    --reference data/HXB2R.fasta \
    --amplicons data/amplicons.tsv
```

**Benefits:**
- Fast local development (no container rebuilds)
- Easy tool version updates (edit YAML)
- Automatic environment caching
- Works on systems without Docker/Singularity

**Files affected:**
- `env/*.yml` (new directory and files)
- `nextflow.config` - add conda profile
- All module files in `modules/local/*.nf` - add `conda` directive

### 5. Dockerfile Improvements

**Goal:** Build and install the local codonyat package instead of pulling a fixed version from PyPI, ensuring the container uses current code.

**Implementation:**

Update `Dockerfile` to copy and install local package:

```dockerfile
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Install system dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bowtie2 \
        default-jre-headless \
        fastp \
        fastqc \
        picard-tools \
        samtools \
    && mkdir -p /usr/picard \
    && ln -sf "$(find /usr/share/java -name 'picard*.jar' | head -n 1)" /usr/picard/picard.jar \
    && python -m pip install --upgrade pip setuptools wheel \
    && rm -rf /var/lib/apt/lists/*

# Copy local package files
COPY pyproject.toml ./
COPY aa_caller/ ./aa_caller/

# Install codonyat from local source + multiqc
RUN pip install --no-cache-dir . multiqc==1.22.2

CMD ["bash"]
```

Create `.dockerignore` to exclude unnecessary files from build context:

```
# Version control
.git/
.github/
.gitignore

# Python
__pycache__/
*.pyc
*.pyo
*.pyd
.pytest_cache/
*.egg-info/

# Build artifacts
dist/
build/

# Nextflow
work/
results/
.nextflow/
.nextflow.log*

# Documentation and tests
docs/
tests/
README.md
LICENSE

# Data (use mounted volumes instead)
data/

# Environment files
env/

# Misc
*.md
!pyproject.toml
```

**Changes:**
- Remove `codonyat==1.0.0` PyPI install
- Copy `pyproject.toml` and `aa_caller/` directory
- Install with `pip install .` (editable install not needed in container)
- Keep multiqc from PyPI since it's not part of this project

**Benefits:**
- Container always uses repository code
- No version mismatch between code and container
- Faster development iteration
- Single source of truth for package version

**Files affected:**
- `Dockerfile`
- `.dockerignore` (new file)

## Testing Strategy

After implementation, validate each improvement:

1. **Retry strategy:**
   - Simulate transient failure (kill process mid-execution)
   - Verify automatic retry occurs
   - Check logs show retry attempts

2. **Input validation:**
   - Test with missing samplesheet → fails fast with clear error
   - Test with invalid CSV format → fails with column error
   - Test with non-existent FASTQ → fails with file-not-found error
   - Test with duplicate sample IDs → fails with duplicate error

3. **Error handling:**
   - Test with invalid reference FASTA → clear error from bowtie2-build
   - Test with wrong protein name → clear error from codonyat
   - Test with empty FASTQ → clear error from process

4. **Conda profile:**
   - Run test data with `-profile conda`
   - Verify all tools execute correctly
   - Check environment caching works

5. **Dockerfile:**
   - Build container: `docker build -t codonyat:latest .`
   - Verify codonyat command available
   - Check version matches local code
   - Run test pipeline with new container

## Implementation Order

Recommended order for minimal interdependence:

1. **Input validation** (independent, quick win)
2. **Retry strategy** (independent, config-only)
3. **Error handling** (independent, per-module)
4. **Dockerfile improvements** (independent)
5. **Conda profile** (independent, but test after all modules updated)

All improvements are independent and can be implemented in parallel or any order.

## Files Summary

**New files:**
- `lib/Utils.groovy` - input validation functions
- `env/fastqc.yml` - FastQC conda environment
- `env/fastp.yml` - fastp conda environment
- `env/bowtie2.yml` - bowtie2 conda environment
- `env/samtools.yml` - samtools conda environment
- `env/picard.yml` - Picard conda environment
- `env/codonyat.yml` - codonyat conda environment
- `env/multiqc.yml` - MultiQC conda environment
- `.dockerignore` - exclude files from Docker build

**Modified files:**
- `conf/base.config` - add retry directives
- `main.nf` - add validation calls
- `nextflow.config` - add conda profile
- `Dockerfile` - install local package
- `modules/local/bowtie2.nf` - add error handling + conda directive
- `modules/local/codonyat.nf` - add error handling + conda directive
- `modules/local/fastp.nf` - add error handling + conda directive
- `modules/local/fastqc.nf` - add conda directive
- `modules/local/samtools.nf` - add error handling + conda directive
- `modules/local/dedup.nf` - add error handling + conda directive
- `modules/local/multiqc.nf` - add conda directive

## Success Criteria

Implementation is complete when:

1. ✅ Transient failures automatically retry up to configured limits
2. ✅ Invalid samplesheets are caught before pipeline starts
3. ✅ Process failures produce clear, actionable error messages
4. ✅ Pipeline runs successfully with `-profile conda`
5. ✅ Docker container installs local package, not PyPI version
6. ✅ All existing tests pass with new improvements
7. ✅ Test pipeline completes successfully with both docker and conda profiles

## Future Considerations

Items explicitly deferred for simplicity:

- **Per-process containers:** Keep monolithic container for easier deployment
- **nf-core standardization:** Overkill for a pipeline this size
- **Resource auto-scaling:** Current fixed resources work for target workloads
- **Publishing to nf-core modules:** Not planning public distribution

These can be revisited if requirements change.
