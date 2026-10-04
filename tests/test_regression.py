"""Byte-for-byte regression against codonyat 1.0.3 outputs (see tests/data/regression_v1_0_3)."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from aa_caller.runner import call_variants

DATA = Path(__file__).parent / "data" / "regression_v1_0_3"
PROTEINS = ["PR", "RT", "INT"]
SAMPLES = ["regression", "edge"]


def _stage(sample: str, dest: Path) -> Path:
    """Copy the SAM fixture into *dest* as ``<sample>.sam`` (decompressing if needed)."""
    target = dest / f"{sample}.sam"
    gz = DATA / f"{sample}.sam.gz"
    if gz.exists():
        with gzip.open(gz, "rb") as src, open(target, "wb") as out:
            shutil.copyfileobj(src, out)
    else:
        shutil.copy(DATA / f"{sample}.sam", target)
    shutil.copy(DATA / "reference.fasta", dest / "reference.fasta")
    shutil.copy(DATA / "amplicons.tsv", dest / "amplicons.tsv")
    return target


def _expected(sample: str, protein: str, ext: str) -> bytes:
    return gzip.decompress((DATA / "expected" / f"{sample}.{protein}.{ext}.gz").read_bytes())


@pytest.mark.parametrize("sample", SAMPLES)
@pytest.mark.parametrize("protein", PROTEINS)
def test_single_protein_outputs_identical_to_v1_0_3(tmp_path, monkeypatch, sample, protein):
    _stage(sample, tmp_path)
    monkeypatch.chdir(tmp_path)
    result = call_variants(f"{sample}.sam", "reference.fasta", "amplicons.tsv", protein=protein)
    assert result.csv_path.read_bytes() == _expected(sample, protein, "tsv")
    assert result.xml_path.read_bytes() == _expected(sample, protein, "xml")


@pytest.mark.parametrize("sample", SAMPLES)
def test_cli_single_protein_identical_to_v1_0_3(tmp_path, sample):
    _stage(sample, tmp_path)
    subprocess.run(
        [sys.executable, "-m", "aa_caller", f"{sample}.sam", "reference.fasta", "amplicons.tsv", "--protein", "RT"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    assert (tmp_path / f"{sample}.sam.tsv").read_bytes() == _expected(sample, "RT", "tsv")
    assert (tmp_path / f"{sample}.sam.xml").read_bytes() == _expected(sample, "RT", "xml")


@pytest.mark.parametrize("sample", SAMPLES)
def test_multi_protein_matches_single_protein_runs(tmp_path, monkeypatch, sample):
    """One pass over PR,RT,INT gives the rows of the three single-protein runs, in order."""
    _stage(sample, tmp_path)
    monkeypatch.chdir(tmp_path)
    result = call_variants(f"{sample}.sam", "reference.fasta", "amplicons.tsv", protein="PR,RT,INT")
    assert result.proteins == PROTEINS

    lines = result.csv_path.read_text().splitlines()
    expected_lines = []
    for i, protein in enumerate(PROTEINS):
        single = _expected(sample, protein, "tsv").decode().splitlines()
        expected_lines.extend(single if i == 0 else single[1:])
    assert lines == expected_lines

    root = ET.parse(result.xml_path).getroot()
    assert [p.get("name") for p in root.findall("Protein")] == PROTEINS
    for protein_elem in root.findall("Protein"):
        single_root = ET.fromstring(_expected(sample, protein_elem.get("name"), "xml"))
        assert [ET.tostring(e) for e in protein_elem] == [ET.tostring(e) for e in single_root]


def test_all_proteins_follow_header_order(tmp_path, monkeypatch):
    _stage("edge", tmp_path)
    monkeypatch.chdir(tmp_path)
    result = call_variants("edge.sam", "reference.fasta", "amplicons.tsv", protein="all")
    assert result.proteins == ["PR", "RT", "INT"]


def test_gzipped_sam_gives_same_counts(tmp_path, monkeypatch):
    sam = _stage("regression", tmp_path)
    gz = tmp_path / "regression.sam.gz"
    with open(sam, "rb") as src, gzip.open(gz, "wb") as out:
        shutil.copyfileobj(src, out)
    monkeypatch.chdir(tmp_path)
    result = call_variants(
        "regression.sam.gz", "reference.fasta", "amplicons.tsv", protein="RT", csv_path="gz.tsv", xml_path="gz.xml"
    )
    got = [line.split("\t", 1)[1] for line in result.csv_path.read_text().splitlines()[1:]]
    want = [line.split("\t", 1)[1] for line in _expected("regression", "RT", "tsv").decode().splitlines()[1:]]
    assert got == want


def test_bam_gives_same_counts(tmp_path, monkeypatch):
    pysam = pytest.importorskip("pysam")
    sam = _stage("regression", tmp_path)
    bam = tmp_path / "regression.bam"
    with pysam.AlignmentFile(str(sam), "r") as src, pysam.AlignmentFile(str(bam), "wb", template=src) as out:
        for rec in src:
            out.write(rec)
    monkeypatch.chdir(tmp_path)
    result = call_variants(
        "regression.bam", "reference.fasta", "amplicons.tsv", protein="RT", csv_path="bam.tsv", xml_path="bam.xml"
    )
    got = [line.split("\t", 1)[1] for line in result.csv_path.read_text().splitlines()[1:]]
    want = [line.split("\t", 1)[1] for line in _expected("regression", "RT", "tsv").decode().splitlines()[1:]]
    assert got == want
