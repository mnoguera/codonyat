# aa-caller

A standalone Python package that exposes the same RT variant calling workflow as `stand_alone_fullkkwr2.py`, including:

- SAM parsing that honors amplicon labels and strand orientation.
- Ratio balancing, entropy tracking, and TSV/XML exporters mirroring the legacy Perl outputs.
- Configurable CLI flags for strand ratio bounds and entropy sensitivity.
 
## Highlights

- Produces TSV/XML diagnostics identical to the legacy Perl script while adding Python objects (`SamContainer`, `FullReference`, etc.) that downstream tooling can import.
- Validates every input file before parsing to surface malformed SAM/FASTA/amplicon data early.
- Logs per-position entropy and strand balance for easier debugging in CI or local runs.

## Installation

```bash
pip install ./aa_caller
```

Or publish the package (e.g., via `twine`/PyPI) and install it like any other dependency.

## CLI usage

Once installed, the `aa-caller` entry point is available:

```bash
aa-caller /path/to/sample.sam /path/to/reference.fasta /path/to/amplicons.tsv
```

Supply `--ratio-upper`, `--ratio-lower`, and `--entropy-threshold` to tune the balancing heuristics.

The CLI writes `[sam-file].tsv` (columns: FILE, REFERENCE, PROTEIN, VARIANT, POSITION, FREQ, FWCOV, RVCOV, TOTALCOV, RATIO) and `[sam-file].xml` (per-position `<Depth>`, `<FwCover>`, `<RvCover>`, `<Variants>`).

## Package API

Import `aa_caller` to reuse the core objects:

```python
from aa_caller import SamContainer, FullReference, parse_amplicons
```

The `SamContainer` constructor still accepts `ratio_upper`, `ratio_lower`, and `entropy_threshold` so you can reuse the balancing logic in scripts.

## Development

Install the repository with the optional dev tooling so your local environment matches CI:

```bash
pip install --upgrade pip
pip install -e .[dev]
```

Now you can run the same checks that land in [.github/workflows/python-tests.yml](.github/workflows/python-tests.yml#L1-L27):

```bash
ruff check .
python -m pytest
```

The workflow installs `pytest` and `ruff`, runs the linter, and then executes the pytest suite on every push/PR against `main`.

## Testing

Run the upstream validation helpers with `pytest`:

```bash
python -m pytest
```

Keep the code tidy with `ruff` before committing:

```bash
ruff check .
```
