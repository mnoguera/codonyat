# Benchmarks

Development tooling only; not part of the installed package.

## Peak memory and run time, 1.0.3 vs 1.1.0

Data: synthetic reads over HXB2 pol from `make_synthetic_sam.py` (~250 nt reads,
substitutions, 1/3-base indels, soft clips, both strands, ~1% unmapped), using the
HXB2 reference with `PR(Protease):2253-2549;RT(Reverse Transcriptase):2550-3869;INT(Integrase):4230-5093`.

```bash
python benchmarks/make_synthetic_sam.py reference.fasta syn_1000000.sam --reads 1000000
/usr/bin/time -v codonyat syn_1000000.sam reference.fasta amplicons.tsv --protein RT
```

Measured on 4 Oct 2026 (Debian 13, Python 3.13.5, one core), "Maximum resident set size" from GNU time:

| Alignments | Protein(s) | 1.0.3 peak RSS | 1.0.3 time | 1.1.0 peak RSS | 1.1.0 time | Outputs |
|---:|---|---:|---:|---:|---:|---|
| 20,000 | RT | 0.48 GB | 3.7 s | n/m | n/m | |
| 50,000 | RT | 1.13 GB | 14.7 s | 64 MB | 1.2 s | identical |
| 100,000 | RT | 2.20 GB | 30.9 s | 69 MB | 2.1 s | identical (PR, INT also identical) |
| 200,000 | RT | 4.34 GB | 60.1 s | n/m | n/m | identical |
| 1,000,000 | RT | ~22 GB (extrapolated, ~22 KB/alignment) | ~5 min (extrapolated) | **83 MB** | 17.2 s | |
| 1,000,000 | PR+RT+INT (`all`), one pass | ~22 GB × 3 runs | | **112 MB** | 28.7 s | RT rows identical to the single-protein run |

1.0.3 was not run at 1M alignments because the box has 15 GB of RAM; its memory
grows linearly (0.48 → 1.13 → 2.20 → 4.34 GB), so the 1M figure is extrapolated.
"identical" means `cmp` found the TSV and XML byte-identical to 1.0.3.
