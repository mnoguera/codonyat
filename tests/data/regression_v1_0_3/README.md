# Regression fixtures (codonyat 1.0.3)

Expected outputs in `expected/` were produced by **codonyat 1.0.3** (the
load-everything implementation, upstream `main` @ `bac40d5`) with:

```
codonyat <name>.sam reference.fasta amplicons.tsv --protein <PR|RT|INT>
```

run from the directory containing `<name>.sam` (so the `FILE` column and the
XML `sample` attribute are just `<name>.sam`).

- `regression.sam.gz`: 2,000 synthetic reads over HXB2 pol from
  `benchmarks/make_synthetic_sam.py reference.fasta regression.sam --reads 2000 --read-length 150 --seed 7`
  (substitutions, 1/3-base indels, soft clips, both strands, unmapped records,
  `Amp_<n>` tags, `_<n>` multiplicity suffixes).
- `edge.sam`: hand-written records for CIGAR edge cases (D, I, N, H, P, `=`/`X`,
  mid-read soft clip, off-frame starts, reads crossing protein boundaries,
  too-short and unmapped records, multiplicity suffixes).

`tests/test_regression.py` checks that the streaming implementation
reproduces these files byte for byte.
