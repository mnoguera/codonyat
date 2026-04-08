# codonyat

codon-yat — A codon-aware amino acid variant typer from SAM alignments of viral NGS data.

A standalone Python package and Nextflow pipeline for amino acid variant calling from viral amplicon sequencing, including:

- SAM parsing that honors amplicon labels and strand orientation.
- CIGAR-aware codon extraction that correctly handles insertions, deletions, and soft clips.
- Ratio balancing, entropy tracking, and TSV/XML exporters mirroring the legacy Perl outputs.
- Configurable CLI flags for strand ratio bounds, entropy sensitivity, and target protein.
- A complete Nextflow DSL2 pipeline: QC, trimming, decontamination, alignment, optional dedup, and variant calling.

## Highlights

- Produces consistent TSV/XML diagnostics while exposing Python objects (`SamContainer`, `FullReference`, etc.) that downstream tooling can import.
- Validates every input file before parsing to surface malformed SAM/FASTA/amplicon data early.
- Logs per-position entropy and strand balance for easier debugging in CI or local runs.
- Supports any protein target via `--protein` (defaults to RT).

## Installation

```bash
pip install .
```

Or with Docker:

```bash
docker build -t codonyat:latest .
```

## CLI usage

Once installed, the `codonyat` entry point is available:

```bash
codonyat sample.sam reference.fasta amplicons.tsv
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--protein` | `RT` | Target protein name (must match a FASTA header annotation) |
| `--ratio-upper` | `3.162` | Upper bound for strand-ratio balancing |
| `--ratio-lower` | `0.316` | Lower bound for strand-ratio balancing |
| `--entropy-threshold` | `0.0` | Minimum Shannon entropy to mark a position as diverse |

### Output

The CLI writes two files alongside the input SAM:

- `[sam-file].tsv` — columns: FILE, REFERENCE, PROTEIN, VARIANT, POSITION, FREQ, FWCOV, RVCOV, TOTALCOV, RATIO
- `[sam-file].xml` — per-position `<Depth>`, `<FwCover>`, `<RvCover>`, `<Variants>`

### Reference FASTA format

The reference FASTA header must include protein annotations in the format `Name(Description):start-end`, separated by semicolons:

```
>K03455|HIVHXB2CG PR(Protease):2253-2549;RT(Reverse Transcriptase):2550-3869;INT(Integrase):4230-5093
TGGAAGGGCTAATTCACTCCCAACGAAGAC...
```

## Nextflow pipeline

The repository includes a Nextflow DSL2 pipeline that wraps codonyat with upstream processing steps.

### Pipeline steps

```
FASTQ → FastQC (raw) → fastp (trim/filter) → FastQC (trimmed)
  → bowtie2 (align, --very-sensitive-local --no-unal)
  → samtools sort → [Picard MarkDuplicates] → BAM → SAM → codonyat
  → MultiQC
```

| Step | Tool | Purpose |
|------|------|---------|
| Quality control | FastQC | Per-read QC metrics on raw and trimmed reads |
| Trimming | fastp | Adapter removal, quality/length filtering |
| Decontamination | bowtie2 `--no-unal` | Drops reads that don't align to the reference |
| Alignment | bowtie2 `--very-sensitive-local` | Sensitive local alignment for diverse viral populations |
| Sorting | samtools sort | Coordinate-sorted BAM |
| Deduplication | Picard MarkDuplicates | Optional (off by default for amplicon data) |
| Variant calling | codonyat | Codon-aware amino acid variant typing |
| Report | MultiQC | Aggregated QC report |

### Quick start

```bash
# Build the codonyat container
docker build -t codonyat:latest .

# Run with your data
nextflow run . -profile docker \
    --input samplesheet.csv \
    --reference data/HXB2R.fasta \
    --amplicons data/amplicons.tsv

# Run with bundled test data
nextflow run . -profile docker,test
```

### Samplesheet format

A CSV file with columns `sample_id`, `fastq_1`, and `fastq_2` (leave `fastq_2` empty for single-end):

```csv
sample_id,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### Pipeline parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Samplesheet CSV |
| `--reference` | (required) | FASTA reference with protein annotations in the header |
| `--amplicons` | (required) | Amplicon definitions TSV |
| `--protein` | `RT` | Target protein for variant calling |
| `--outdir` | `./results` | Output directory |
| `--skip_dedup` | `true` | Skip deduplication (default for amplicon data) |
| `--skip_fastqc` | `false` | Skip FastQC steps |
| `--bowtie2_args` | `--very-sensitive-local --no-unal` | bowtie2 alignment flags |
| `--fastp_args` | `--qualified_quality_phred 20 --length_required 50 --detect_adapter_for_pe` | fastp trimming flags |
| `--ratio_upper` | `3.162` | Upper strand-ratio bound (passed to codonyat) |
| `--ratio_lower` | `0.316` | Lower strand-ratio bound (passed to codonyat) |
| `--entropy_threshold` | `0.0` | Minimum Shannon entropy (passed to codonyat) |

### Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `test` | Use bundled test data from `data/` |

### Output structure

```
results/
├── fastqc_raw/           # FastQC reports on raw reads
├── fastp/                # Trimming reports and logs
├── fastqc_trimmed/       # FastQC reports on trimmed reads
├── bowtie2_align/        # Alignment BAMs and logs
├── samtools_sort/        # Sorted BAMs
├── codonyat/             # Variant TSV and XML per sample
└── multiqc/              # Aggregated QC report
```

## Package API

Import `aa_caller` to reuse the core objects:

```python
from aa_caller import SamContainer, FullReference, parse_amplicons
```

### Python wrapper

Use `call_variants` to run the full pipeline from Python without touching the CLI:

```python
from pathlib import Path
from aa_caller import call_variants

result = call_variants(
    sam_path=Path("reads.sam"),
    reference_path=Path("reference.fasta"),
    amplicons_path=Path("amps.tsv"),
    protein="RT",
    ratio_upper=3.2,
)

print(result.csv_path, result.xml_path)
print(result.container.variants.keys())
```

If you already have an `argparse.Namespace` or mapping of the CLI arguments, `call_variants_from_args` adapts them directly:

```python
from argparse import Namespace

args = Namespace(
    sam_file="reads.sam",
    reference_file="reference.fasta",
    amplicons_file="amps.tsv",
    protein="RT",
    ratio_upper=3.2,
)

call_variants_from_args(args)
```

### CLI wrapper

A lightweight runner at `codonyat-runner` exposes the same arguments as `call_variants_from_args`:

```bash
codonyat-runner reads.sam reference.fasta amps.tsv --protein RT --ratio-upper 3.1 --csv-path results.tsv
```

## Development

Install the repository with the optional dev tooling so your local environment matches CI:

```bash
pip install --upgrade pip
pip install -e .[dev]
```

Now you can run the same checks that land in [.github/workflows/python-tests.yml](.github/workflows/python-tests.yml):

```bash
ruff check .
python -m pytest
```

## Project structure

```
codonyat/
├── aa_caller/            # Python package
│   ├── cli.py            # CLI entry point
│   ├── runner.py         # Python API (call_variants, call_variants_from_args)
│   ├── container.py      # SamContainer — variant aggregation engine
│   ├── sam.py            # SamEntry — CIGAR-aware SAM record parser
│   ├── reference.py      # FullReference — FASTA + protein annotation loader
│   ├── models.py         # Protein, Amplicon, Variant, Qual dataclasses
│   ├── validators.py     # Input file validation
│   ├── genetic_code.py   # Codon translation table
│   └── constants.py      # Default thresholds
├── modules/local/        # Nextflow process modules
├── conf/                 # Nextflow resource configuration
├── data/                 # Test data (HXB2R reference, amplicons, sample SAM)
├── tests/                # pytest suite (unit + integration)
├── main.nf               # Nextflow pipeline entry point
├── nextflow.config       # Pipeline parameters and profiles
├── Dockerfile            # Container build for codonyat
└── pyproject.toml        # Python package metadata
```

## License

MIT
