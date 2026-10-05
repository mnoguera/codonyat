# Changelog

## 1.1.0 (unreleased)

### Changed
- **Streaming input.** Alignment records are processed one at a time instead of
  being loaded into memory first. Peak memory no longer grows with the number of
  alignments: on 1,000,000 synthetic HXB2 pol alignments peak RSS is **~83 MB**
  for RT (1.0.3: ~22 GB extrapolated from 4.3 GB at 200,000 alignments, about
  22 KB per alignment) and **~112 MB** for PR+RT+INT in one pass. Run time drops
  from O(positions × reads) to O(reads): 17 s for RT at 1M alignments
  (1.0.3: ~5 min extrapolated). See `benchmarks/README.md`.
- **Single-protein outputs are byte-for-byte identical to 1.0.3**, checked by
  `tests/test_regression.py` against outputs produced by 1.0.3.
- `call_variants()` / the CLIs no longer call `SamContainer.load()`, so
  `result.container.reads` is now empty. The new `container.records_seen` and
  `container.mapped_records` give the counts. `load()` still works for code that
  needs the reads in memory, and `calculate_variant_frequencies()` uses them if
  they were loaded.

### Added
- **Multi-protein runs in one pass:** `--protein PR,RT,INT` or `--protein all`
  (every protein annotated in the reference header, in header order). The TSV
  has one header and the rows of each protein in the requested order; the XML
  wraps each protein's positions in a `<Protein name="...">` element. A single
  protein produces exactly the same TSV/XML layout as before.
- Python API: `SamContainer.calculate_variant_frequencies_multi()`,
  `write_csv_proteins()`, `write_xml_proteins()`, `container.protein_variants`,
  `resolve_proteins()`, `parse_protein_spec()`, `iter_sam_lines()`;
  `call_variants(protein=...)` accepts a comma-separated string, a list or `"all"`,
  and `VariantCallResult.proteins` lists the proteins analysed.
- Input can be SAM, gzip-compressed SAM (`.sam.gz`) or BAM (`.bam`, needs the
  optional extra: `pip install 'codonyat[bam]'`).
- A record whose SEQ is shorter than its CIGAR now raises a clear `ValueError`
  (1.0.3 failed with an `IndexError`).
- `benchmarks/make_synthetic_sam.py` (deterministic synthetic SAM generator) and
  `benchmarks/README.md` (how the memory numbers were measured).

## 1.0.3
- Python-package-only repository; release-based PyPI publishing.
