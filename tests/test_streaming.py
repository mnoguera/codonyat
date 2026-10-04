"""Streaming behaviour: flat memory, legacy load() path, protein selection."""

from __future__ import annotations

import random
import tracemalloc
from pathlib import Path

import pytest

from aa_caller.constants import DEFAULT_ENTROPY_THRESHOLD, DEFAULT_RATIO_LOWER, DEFAULT_RATIO_UPPER
from aa_caller.container import SamContainer
from aa_caller.reference import FullReference, parse_protein_spec, resolve_proteins

DATA = Path(__file__).parent / "data" / "regression_v1_0_3"


def _synthetic_sam(path: Path, n_reads: int, seed: int = 3) -> Path:
    """Simple 150M reads over HXB2 PR/RT, enough to exercise the streaming path.

    Variation is limited to one fixed site so codon diversity (which, unlike
    the read count, legitimately sizes the result tables) stays constant.
    """
    ref = "".join(line.strip() for line in (DATA / "reference.fasta").read_text().splitlines()[1:])
    rng = random.Random(seed)
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\n")
        for i in range(n_reads):
            pos = rng.randint(2200, 3700)
            seq = list(ref[pos - 1 : pos + 149])
            if pos <= 2700 < pos + 150 and i % 3 == 0:
                seq[2700 - pos] = "N"
            seq = "".join(seq)
            flag = 16 if i % 2 else 0
            fh.write(f"r{i}:Amp_1\t{flag}\tref\t{pos}\t60\t150M\t*\t0\t0\t{seq}\t{'I' * 150}\n")
    return path


def _peak_bytes(sam: Path, reference: FullReference) -> int:
    container = SamContainer(sam, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
    tracemalloc.start()
    container.calculate_variant_frequencies_multi(["PR", "RT"], reference)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert container.mapped_records > 0
    return peak


def test_memory_does_not_grow_with_alignment_count(tmp_path):
    reference = FullReference(DATA / "reference.fasta")
    small = _peak_bytes(_synthetic_sam(tmp_path / "small.sam", 2_000), reference)
    large = _peak_bytes(_synthetic_sam(tmp_path / "large.sam", 20_000), reference)
    # 10x more alignments must not mean ~10x more memory (1.0.3 grew linearly, ~24 KB/alignment,
    # i.e. ~480 MB for 20,000 alignments).
    assert large < small * 1.5 + 1_000_000
    assert large < 20_000_000


def test_legacy_load_path_gives_same_results(tmp_path):
    reference = FullReference(DATA / "reference.fasta")
    sam = _synthetic_sam(tmp_path / "reads.sam", 500)
    streamed = SamContainer(sam, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
    streamed.calculate_variant_frequencies("RT", reference)
    loaded = SamContainer(sam, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
    loaded.load()
    assert len(loaded.reads) == 500
    loaded.calculate_variant_frequencies("RT", reference)
    streamed.write_csv("x", reference, tmp_path / "a.tsv", protein_name="RT")
    loaded.write_csv("x", reference, tmp_path / "b.tsv", protein_name="RT")
    assert (tmp_path / "a.tsv").read_bytes() == (tmp_path / "b.tsv").read_bytes()


def test_parse_protein_spec():
    assert parse_protein_spec("RT") == ["RT"]
    assert parse_protein_spec(" PR, RT ,PR,INT ") == ["PR", "RT", "INT"]
    assert parse_protein_spec(["RT", "PR"]) == ["RT", "PR"]
    assert parse_protein_spec("all") == ["all"]
    with pytest.raises(ValueError):
        parse_protein_spec(" , ")
    with pytest.raises(ValueError):
        parse_protein_spec("all,RT")


def test_resolve_proteins():
    reference = FullReference(DATA / "reference.fasta")
    assert resolve_proteins("ALL", reference) == ["PR", "RT", "INT"]
    assert resolve_proteins("INT,PR", reference) == ["INT", "PR"]
    with pytest.raises(ValueError, match="GAG"):
        resolve_proteins("RT,GAG", reference)


def test_multi_requires_known_protein(tmp_path):
    reference = FullReference(DATA / "reference.fasta")
    sam = _synthetic_sam(tmp_path / "reads.sam", 10)
    container = SamContainer(sam, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
    with pytest.raises(ValueError):
        container.calculate_variant_frequencies_multi(["NOPE"], reference)
    with pytest.raises(ValueError):
        container.calculate_variant_frequencies_multi([], reference)


def test_seq_shorter_than_cigar_is_reported(tmp_path):
    reference = FullReference(DATA / "reference.fasta")
    sam = tmp_path / "bad.sam"
    sam.write_text("r1\t0\tref\t2253\t60\t30M\t*\t0\t0\tACGT\tIIII\n")
    container = SamContainer(sam, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
    with pytest.raises(ValueError, match="shorter than its CIGAR"):
        container.calculate_variant_frequencies("PR", reference)
    sam.write_text("r1\t0\tref\t2253\t60\t2M1D28M\t*\t0\t0\tACGT\tIIII\n")
    with pytest.raises(ValueError, match="shorter than its CIGAR"):
        container.calculate_variant_frequencies("PR", reference)


def test_cli_multi_protein_and_bad_protein(tmp_path, monkeypatch):
    import shutil
    import sys

    from aa_caller import cli

    sam = _synthetic_sam(tmp_path / "reads.sam", 50)
    shutil.copy(DATA / "reference.fasta", tmp_path / "reference.fasta")
    shutil.copy(DATA / "amplicons.tsv", tmp_path / "amplicons.tsv")
    argv = ["codonyat", str(sam), str(tmp_path / "reference.fasta"), str(tmp_path / "amplicons.tsv")]
    monkeypatch.setattr(sys, "argv", argv + ["--protein", "all"])
    cli.main()
    text = (tmp_path / "reads.sam.tsv").read_text()
    assert "\tPR\t" in text and "\tRT\t" in text
    assert "<Protein name=\"INT\"" in (tmp_path / "reads.sam.xml").read_text()
    monkeypatch.setattr(sys, "argv", argv + ["--protein", "all,RT"])
    with pytest.raises(SystemExit):
        cli.main()
