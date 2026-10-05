# codonyat

codon-yat — A codon-aware amino acid variant typer from SAM alignments of viral NGS data.

A Python package for codon-level amino acid variant calling from SAM alignments of viral amplicon sequencing, including:

- SAM parsing that honors amplicon labels and strand orientation.
- CIGAR-aware codon extraction that correctly handles insertions, deletions, and soft clips.
- Ratio balancing, entropy tracking, and TSV/XML exporters mirroring the legacy Perl outputs.
- Configurable CLI flags for strand ratio bounds, entropy sensitivity, and target protein.

## Highlights

- Produces consistent TSV/XML diagnostics while exposing Python objects (`SamContainer`, `FullReference`, etc.) that downstream tooling can import.
- Validates every input file before parsing to surface malformed SAM/FASTA/amplicon data early.
- Logs per-position entropy and strand balance for easier debugging in CI or local runs.
- Supports any protein target via `--protein` (defaults to RT).

## Installation

```bash
pip install codonyat
```

Requires Python 3.10 or newer. The only runtime dependency is Biopython.

## Running it in a workflow

This repository contains only the Python package. To run codonyat on many samples
(from FASTQ or SAM, with read QC, trimming, alignment and a merged report), use the
separate Nextflow pipeline **[nf-codonyat](https://github.com/sjoclaudi/nf-codonyat)**,
which installs this package from PyPI.

## CLI usage

Once installed, the `codonyat` entry point is available:

```bash
codonyat sample.sam reference.fasta amplicons.tsv
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--protein` | `RT` | Protein(s) to analyse: one name, a comma-separated list (e.g. `PR,RT,INT`, counted in a single pass) or `all` (every protein annotated in the reference header). Names must match the FASTA header annotations |
| `--ratio-upper` | `3.162` | Upper bound for strand-ratio balancing |
| `--ratio-lower` | `0.316` | Lower bound for strand-ratio balancing |
| `--entropy-threshold` | `0.0` | Minimum Shannon entropy to mark a position as diverse |

### Output

The CLI writes two files alongside the input SAM:

- `[sam-file].tsv` — columns: FILE, REFERENCE, PROTEIN, VARIANT, POSITION, FREQ, FWCOV, RVCOV, TOTALCOV, RATIO
- `[sam-file].xml` — per-position `<Depth>`, `<FwCover>`, `<RvCover>`, `<Variants>`

With several proteins, the TSV has one header followed by each protein's rows (in the order requested), and the XML wraps each protein's positions in a `<Protein name="...">` element. Single-protein output is unchanged.

### Input and memory

The input can be SAM, gzip-compressed SAM (`.sam.gz`) or BAM (`.bam`; install the optional extra with `pip install 'codonyat[bam]'`). Records are streamed, so peak memory does not depend on the number of alignments: about 80–110 MB for 1 million alignments (see [benchmarks/README.md](benchmarks/README.md)).

### Reference FASTA format

The reference FASTA header must include protein annotations in the format `Name(Description):start-end`, separated by semicolons:

```
>K03455|HIVHXB2CG PR(Protease):2253-2549;RT(Reverse Transcriptase):2550-3869;INT(Integrase):4230-5093
TGGAAGGGCTAATTCACTCCCAACGAAGAC...
```

## Package API

Import `aa_caller` to reuse the core objects:

```python
from aa_caller import SamContainer, FullReference, parse_amplicons
```

### Python wrapper

Use `call_variants` to run variant calling from Python without touching the CLI:

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

Now you can run the same checks as [.github/workflows/python-tests.yml](.github/workflows/python-tests.yml):

```bash
ruff check .
python -m pytest
python -m build && twine check dist/*
```

## Project structure

```
codonyat/
├── aa_caller/            # Python package
│   ├── cli.py            # CLI entry point (codonyat)
│   ├── runner.py         # Python API (call_variants, call_variants_from_args) and codonyat-runner
│   ├── container.py      # SamContainer — variant aggregation engine
│   ├── sam.py            # SamEntry — CIGAR-aware SAM record parser
│   ├── reference.py      # FullReference — FASTA + protein annotation loader
│   ├── models.py         # Protein, Amplicon, Variant, Qual dataclasses
│   ├── validators.py     # Input file validation
│   ├── genetic_code.py   # Codon translation table
│   └── constants.py      # Default thresholds
├── tests/                # pytest suite (unit, integration, 1.0.3 regression fixtures)
├── benchmarks/           # synthetic SAM generator and memory benchmark notes (dev only)
├── CHANGELOG.md
├── docs/index.html       # Project page (GitHub Pages)
└── pyproject.toml        # Package metadata
```

## Releasing

Publishing to PyPI is done by `.github/workflows/publish-pypi.yml` when a GitHub
release is published:

1. Set `version` in `pyproject.toml` (e.g. `1.0.3`) and merge it to `main`.
2. Create a tag `v<version>` and a GitHub release from it.
3. The workflow builds the sdist and wheel, runs `twine check`, checks that the tag
   matches the package version, and publishes via PyPI trusted publishing
   (GitHub environment `pypi`).

## License

MIT
