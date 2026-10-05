from __future__ import annotations

import csv
import logging
import math
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Dict, List

from .genetic_code import codon_to_aminoacid
from .models import Variant
from .reference import FullReference
from .sam import SamEntry, iter_sam_lines

logger = logging.getLogger(__name__)

# Same patterns as SamEntry, compiled once.
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_AMP_RE = re.compile(r"(Amp_[0-9]+)")

#: How often (in alignment records) streaming progress is logged.
PROGRESS_EVERY = 1_000_000


class _ProteinAccumulator:
    """Streaming codon counter for one protein.

    Counts are keyed by ``(position, codon, amplicon, is_reverse)``. Python
    dicts keep insertion order, so iterating the keys reproduces the order in
    which codons and amplicons were first seen at each position while walking
    the reads in file order. That is exactly the order the original
    load-everything implementation produced, which keeps outputs identical.
    """

    __slots__ = ("name", "start", "end", "counts")

    def __init__(self, name: str, start: int, end: int) -> None:
        self.name = name
        self.start = start
        # ``range(start, end, 3)`` semantics: ``end`` itself is excluded.
        self.end = end
        self.counts: Dict[tuple, int] = {}


def _accumulate_line(line: str, accumulators: Sequence[_ProteinAccumulator]) -> bool:
    """Parse one SAM text record and add its codons to every accumulator.

    Returns False when the record is unmapped (flag 0x4) and was skipped.
    Semantics match :class:`SamEntry` exactly: a read contributes a codon at a
    protein position ``pos`` when ``pos >= POS`` and ``pos + 2`` is still inside
    the reference span covered by its CIGAR; deleted/skipped reference bases
    are reported as ``-``.
    """
    fields = line.strip().split("\t", 11)
    flag = int(fields[1])
    if flag & 4:
        return False
    identifier = fields[0]
    coordinate = int(fields[3])
    cigar = fields[5]
    sequence = fields[9]

    tail = identifier.split("_")[-1]
    occurrences = max(1, int(tail)) if tail.isdigit() else 1
    amp_match = _AMP_RE.search(identifier)
    amplicon = amp_match.group(1) if amp_match else "Amp_NONE"
    reverse = bool(flag & 16)

    # Walk the CIGAR once. ``refmap`` (ref offset -> read index, -1 = gap) is
    # only built when the alignment is not a simple linear block.
    ops = _CIGAR_RE.findall(cigar)
    reflen = 0
    lead = 0
    linear = True
    seen_ref = False
    clipped_after_ref = False
    for length_str, op in ops:
        length = int(length_str)
        if op in "M=X":
            if clipped_after_ref:
                linear = False
            seen_ref = True
            reflen += length
        elif op in "DN":
            linear = False
            seen_ref = True
            reflen += length
        elif op in "IS":
            if seen_ref:
                clipped_after_ref = True
            else:
                lead += length
    if reflen < 3:
        return True
    last_codon_start = coordinate + reflen - 3  # pos + 2 must stay in the span

    refmap: List[int] | None = None
    if not linear:
        refmap = []
        read_pos = 0
        for length_str, op in ops:
            length = int(length_str)
            if op in "M=X":
                refmap.extend(range(read_pos, read_pos + length))
                read_pos += length
            elif op in "IS":
                read_pos += length
            elif op in "DN":
                refmap.extend([-1] * length)

    for acc in accumulators:
        start = acc.start
        if coordinate > start:
            first = start + -(-(coordinate - start) // 3) * 3
        else:
            first = start
        stop = min(acc.end, last_codon_start + 1)
        if first >= stop:
            continue
        counts = acc.counts
        if refmap is None:
            offset = lead - coordinate
            for pos in range(first, stop, 3):
                i = pos + offset
                codon = sequence[i : i + 3]
                if len(codon) < 3:
                    raise ValueError(
                        f"SAM record {identifier!r}: SEQ is shorter than its CIGAR"
                    )
                key = (pos, codon, amplicon, reverse)
                counts[key] = counts.get(key, 0) + occurrences
        else:
            for pos in range(first, stop, 3):
                o = pos - coordinate
                a, b, c = refmap[o], refmap[o + 1], refmap[o + 2]
                try:
                    codon = (
                        ("-" if a < 0 else sequence[a])
                        + ("-" if b < 0 else sequence[b])
                        + ("-" if c < 0 else sequence[c])
                    )
                except IndexError as exc:
                    raise ValueError(
                        f"SAM record {identifier!r}: SEQ is shorter than its CIGAR"
                    ) from exc
                key = (pos, codon, amplicon, reverse)
                counts[key] = counts.get(key, 0) + occurrences
    return True


class SamContainer:
    """Aggregates SAM records (streamed) and computes codon variant statistics per protein."""

    def __init__(
        self,
        sam_path: Path,
        ratio_upper: float,
        ratio_lower: float,
        entropy_threshold: float,
    ) -> None:
        """Initialize the container with the SAM path, ratio window, and entropy threshold."""
        self.sam_path = sam_path
        self.reads: List[SamEntry] = []
        self.variants: Dict[int, dict] = {}
        self.protein_variants: Dict[str, Dict[int, dict]] = {}
        self.records_seen = 0
        self.mapped_records = 0
        self.ratio_upper = ratio_upper
        self.ratio_lower = ratio_lower
        self.entropy_threshold = entropy_threshold

    def load(self) -> None:
        """Load every mapped read into ``self.reads`` (legacy, memory-hungry).

        Kept for API compatibility. The CLI and :func:`call_variants` no longer
        call it: :meth:`calculate_variant_frequencies` streams the SAM file
        directly when nothing has been loaded, which keeps memory flat
        regardless of the number of alignments.
        """
        for line in iter_sam_lines(self.sam_path):
            entry = SamEntry(line)
            if entry.is_mapped():
                self.reads.append(entry)

    def _record_lines(self) -> Iterable[str]:
        """Records to aggregate: the loaded reads if :meth:`load` was used, else a stream."""
        if self.reads:
            return (read.raw for read in self.reads)
        return iter_sam_lines(self.sam_path)

    def calculate_variant_frequencies(self, protein_name: str, reference: FullReference) -> None:
        """Walk every codon position of the requested protein and build variant histories.

        Equivalent to ``calculate_variant_frequencies_multi([protein_name], reference)``.
        """
        self.calculate_variant_frequencies_multi([protein_name], reference)

    def calculate_variant_frequencies_multi(
        self, protein_names: Sequence[str], reference: FullReference
    ) -> None:
        """Count codons for several proteins in a single streaming pass over the SAM file.

        Results are stored per protein in ``self.protein_variants``;
        ``self.variants`` points at the first requested protein (the only one in
        single-protein runs), as before.
        """
        if not protein_names:
            raise ValueError("At least one protein is required")
        accumulators: List[_ProteinAccumulator] = []
        for name in protein_names:
            if name not in reference.proteins:
                raise ValueError(f"Protein {name} not found in reference")
            protein = reference.proteins[name]
            accumulators.append(
                _ProteinAccumulator(name, protein.start_coordinate, protein.end_coordinate)
            )

        records = 0
        mapped = 0
        for line in self._record_lines():
            records += 1
            if _accumulate_line(line, accumulators):
                mapped += 1
            if records % PROGRESS_EVERY == 0:
                logger.info("Processed %s alignment records", f"{records:,}")
        self.records_seen = records
        self.mapped_records = mapped
        logger.info("Streamed %s records (%s mapped)", f"{records:,}", f"{mapped:,}")

        self.protein_variants = {}
        for acc in accumulators:
            self.variants = self._build_entries(acc)
            for pos in self.variants:
                self.calculate_ratios_nucleotide(pos)
                self.calculate_ratios_aminoacid(pos)
                self.calculate_shannon_entropy_single_position(pos)
                self.log_position_stats(pos)
            self.protein_variants[acc.name] = self.variants
            acc.counts = {}
        self.variants = self.protein_variants[accumulators[0].name]

    @staticmethod
    def _build_entries(acc: _ProteinAccumulator) -> Dict[int, dict]:
        """Turn streamed counts into the per-position structure used by the writers."""
        entries: Dict[int, dict] = {}
        for pos in range(acc.start, acc.end, 3):
            entries[pos] = {
                "variants": {},
                "fw_depth": 0,
                "rv_depth": 0,
                "depth": 0,
                "amplicons": {},
            }
        for (pos, codon, amplicon_label, reverse), n in acc.counts.items():
            entry = entries[pos]
            variant = entry["variants"].get(codon)
            if variant is None:
                variant = entry["variants"][codon] = Variant(codon=codon)
            variant.count += n
            amp_variant = variant.amplicons.setdefault(amplicon_label, {"fw": 0, "rv": 0})
            amp_position = entry["amplicons"].setdefault(
                amplicon_label, {"fw_depth": 0, "rv_depth": 0, "depth": 0}
            )
            if reverse:
                entry["rv_depth"] += n
                variant.rv_reads += n
                amp_variant["rv"] += n
                amp_position["rv_depth"] += n
            else:
                entry["fw_depth"] += n
                variant.fw_reads += n
                amp_variant["fw"] += n
                amp_position["fw_depth"] += n
            amp_position["depth"] += n
            entry["depth"] += n
        return entries

    def calculate_ratios_nucleotide(self, position: int) -> None:
        """Determine strand-specific ratios per variant and flag the balanced reads."""
        entry = self.variants.get(position)
        if not entry or entry["fw_depth"] == 0 or entry["rv_depth"] == 0:
            return
        for variant in entry["variants"].values():
            # compute the per-strand frequency for each codon then derive the ratio
            variant.fw_freq = variant.fw_reads / entry["fw_depth"] if entry["fw_depth"] else 0
            variant.rv_freq = variant.rv_reads / entry["rv_depth"] if entry["rv_depth"] else 0
            variant.ratio = (variant.fw_freq / variant.rv_freq) if variant.rv_freq else 0
            # only mark variants within the balanced ratio window
            if self.ratio_lower <= variant.ratio <= self.ratio_upper:
                variant.balanced_fw += variant.fw_reads
                variant.balanced_rv += variant.rv_reads

    def has_balanced_variants(self, position: int) -> bool:
        """Return True when at least one variant falls within the configured strand-ratio window."""
        entry = self.variants.get(position)
        if not entry:
            return False
        return any(self.ratio_lower <= variant.ratio <= self.ratio_upper for variant in entry["variants"].values())

    def calculate_ratios_aminoacid(self, position: int) -> None:
        """Aggregate amino-acid frequencies across codon variants and capture balance metadata."""
        entry = self.variants.get(position)
        if not entry:
            return
        aa_stats: Dict[str, Dict[str, float]] = {}
        for variant in entry["variants"].values():
            aa = codon_to_aminoacid(variant.codon)
            stats = aa_stats.setdefault(aa, {"count": 0.0, "fw": 0.0, "rv": 0.0})
            stats["count"] += variant.count
            stats["fw"] += variant.fw_reads
            stats["rv"] += variant.rv_reads
        for aa, stats in aa_stats.items():
            fw = stats["fw"]
            rv = stats["rv"]
            ratio = (fw / rv) if rv else 0
            # flag balanced amino-acid counts using the same ratio window as nucleotides
            stats["ratio"] = ratio
            stats["balanced"] = (fw + rv) if self.ratio_lower <= ratio <= self.ratio_upper else 0
        entry["aminoacid_stats"] = aa_stats

    def calculate_shannon_entropy_single_position(self, position: int) -> float:
        """Compute the Shannon entropy for the observed codon mixture at a position."""
        entry = self.variants.get(position)
        if not entry:
            return 0.0
        depth = entry.get("depth", 0)
        if depth == 0:
            entry["shannon_entropy"] = 0.0
            return 0.0
        entropy = 0.0
        for variant in entry["variants"].values():
            if variant.count == 0:
                continue
            p = variant.count / depth
            # standard Shannon formula: sum of -p log(p) for each variant
            entropy -= p * math.log(p)
        entry["shannon_entropy"] = entropy
        entry["entropy_pass"] = entropy >= self.entropy_threshold
        return entropy

    def log_position_stats(self, position: int) -> None:
        """Log strand balance and entropy so downstream debugging is easier."""
        entry = self.variants.get(position)
        if not entry:
            return
        balanced = [v.codon for v in entry["variants"].values() if self.ratio_lower <= v.ratio <= self.ratio_upper]
        entropy = entry.get("shannon_entropy", 0.0)
        logger.debug(
            "Position %s -- depth=%s, fw=%s, rv=%s, entropy=%.3f, balanced=%s",
            position,
            entry["depth"],
            entry["fw_depth"],
            entry["rv_depth"],
            entropy,
            balanced,
        )

    _CSV_FIELDS = [
        "FILE",
        "REFERENCE",
        "PROTEIN",
        "VARIANT",
        "POSITION",
        "FREQ",
        "FWCOV",
        "RVCOV",
        "TOTALCOV",
        "RATIO",
    ]

    def _entries_for(self, protein_name: str) -> Dict[int, dict]:
        """Per-position results for a protein (falls back to ``self.variants``)."""
        return self.protein_variants.get(protein_name, self.variants)

    def _write_csv_rows(self, writer: csv.DictWriter, sample_name: str, reference: FullReference, protein_name: str) -> None:
        for pos, entry in sorted(self._entries_for(protein_name).items()):
            depth = entry["depth"]
            ratio = entry["fw_depth"] / entry["rv_depth"] if entry["rv_depth"] else 0
            for variant in entry["variants"].values():
                # emit each variant using the cached depths and ratios
                writer.writerow(
                    {
                        "FILE": sample_name,
                        "REFERENCE": reference.id,
                        "PROTEIN": protein_name,
                        "VARIANT": variant.codon,
                        "POSITION": pos,
                        "FREQ": round((variant.count / depth * 100) if depth else 0, 3),
                        "FWCOV": entry["fw_depth"],
                        "RVCOV": entry["rv_depth"],
                        "TOTALCOV": depth,
                        "RATIO": round(ratio, 3),
                    }
                )

    def write_csv(self, sample_name: str, reference: FullReference, output_path: Path, *, protein_name: str = "RT") -> None:
        """Emit the variant summary for one protein as a TSV file that mirrors the Perl outputs."""
        self.write_csv_proteins(sample_name, reference, output_path, [protein_name])

    def write_csv_proteins(self, sample_name: str, reference: FullReference, output_path: Path, protein_names: Sequence[str]) -> None:
        """Emit one TSV with a single header and the rows of each protein, in the given order."""
        with open(output_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self._CSV_FIELDS, delimiter="\t")
            writer.writeheader()
            for protein_name in protein_names:
                self._write_csv_rows(writer, sample_name, reference, protein_name)

    @staticmethod
    def _add_positions(parent: ET.Element, entries: Dict[int, dict]) -> None:
        for pos, entry in sorted(entries.items()):
            position_elem = ET.SubElement(parent, "Position", index=str(pos))
            ET.SubElement(position_elem, "Depth").text = str(entry["depth"])
            ET.SubElement(position_elem, "FwCover").text = str(entry["fw_depth"])
            ET.SubElement(position_elem, "RvCover").text = str(entry["rv_depth"])
            variants_elem = ET.SubElement(position_elem, "Variants")
            for variant in entry["variants"].values():
                # describe each variant using XML attributes for diagnostics
                var_elem = ET.SubElement(variants_elem, "Variant", codon=variant.codon)
                var_elem.set("count", str(variant.count))
                var_elem.set("fw_reads", str(variant.fw_reads))
                var_elem.set("rv_reads", str(variant.rv_reads))

    def write_xml(self, sample_name: str, reference: FullReference, output_path: Path) -> None:
        """Dump the per-position state of the (first) protein into XML for diagnostics."""
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        self._add_positions(root, self.variants)
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)

    def write_xml_proteins(self, sample_name: str, reference: FullReference, output_path: Path, protein_names: Sequence[str]) -> None:
        """Multi-protein XML: one ``<Protein name=...>`` element per protein wrapping its positions.

        For a single protein use :meth:`write_xml`, whose layout is unchanged.
        """
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        for protein_name in protein_names:
            protein_elem = ET.SubElement(root, "Protein", name=protein_name)
            self._add_positions(protein_elem, self._entries_for(protein_name))
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)
