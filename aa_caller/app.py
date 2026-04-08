#!/usr/bin/env python3
"""Standalone RT amino acid variant caller translated from stand_alone_fullkkwr2.pl.

This module exposes the full object model replicated from the legacy Perl
implementation so it can be executed within modern Python environments with
identical TSV/XML diagnostics and ratio/entropy helpers that the downstream
pipeline relies on."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


class Qual:
    """Convert SAM quality strings into floating-point scores for filtering."""
    _SOLEXA13 = {
        "@": 0.01,
        "A": 1,
        "B": 2,
        "C": 3,
        "D": 4,
        "E": 5,
        "F": 6,
        "G": 7,
        "H": 8,
        "I": 9,
        "J": 10,
        "K": 11,
        "L": 12,
        "M": 13,
        "N": 14,
        "O": 15,
        "P": 16,
        "Q": 17,
        "R": 18,
        "S": 19,
        "T": 20,
        "U": 21,
        "V": 22,
        "W": 23,
        "X": 24,
        "Y": 25,
        "Z": 26,
        "[": 27,
        "\\": 28,
        "]": 29,
        "^": 30,
        "_": 31,
        "`": 32,
        "a": 33,
        "b": 34,
        "c": 35,
        "d": 36,
        "e": 37,
        "f": 38,
        "g": 39,
        "h": 40,
    }

    _SOLEXA = {
        ";": 0.01,
        "<": 1,
        "=": 2,
        ">": 3,
        "?": 4,
        "@": 5,
        "A": 6,
        "B": 7,
        "C": 8,
        "D": 9,
        "E": 10,
        "F": 11,
        "G": 12,
        "H": 13,
        "I": 14,
        "J": 15,
        "K": 16,
        "L": 17,
        "M": 18,
        "N": 19,
        "O": 20,
        "P": 21,
        "Q": 22,
        "R": 23,
        "S": 24,
        "T": 25,
        "U": 26,
        "V": 27,
        "W": 28,
        "X": 29,
        "Y": 30,
        "Z": 31,
        "[": 32,
        "\\": 33,
        "]": 34,
        "^": 35,
        "_": 36,
        "`": 37,
        "a": 38,
        "b": 39,
        "c": 40,
    }

    _SANGER = {
        "!": 0.01,
        '"': 1,
        "#": 2,
        "$": 3,
        "%": 4,
        "&": 5,
        "'": 6,
        "(": 7,
        ")": 8,
        "*": 9,
        "+": 10,
        ",": 11,
        "-": 12,
        ".": 13,
        "/": 14,
        "0": 15,
        "1": 16,
        "2": 17,
        "3": 18,
        "4": 19,
        "5": 20,
        "6": 21,
        "7": 22,
        "8": 23,
        "9": 24,
        ":": 25,
        ";": 26,
        "<": 27,
        "=": 28,
        ">": 29,
        "?": 30,
        "@": 31,
        "A": 32,
        "B": 33,
        "C": 34,
        "D": 35,
        "E": 36,
        "F": 37,
        "G": 38,
        "H": 39,
        "I": 40,
    }

    def __init__(self, quality_string: str, qual_type: str = "sanger") -> None:
        """Build a numeric representation of the provided FASTQ-style score string."""
        self.quality_string = quality_string
        self.qual_type = qual_type.lower()
        self.quality_array = self._quality_to_numerical_array()

    def _quality_to_numerical_array(self) -> List[float]:
        """Map each character to the configured score table, defaulting to Sanger."""
        mapping = {
            "solexa13": self._SOLEXA13,
            "solexa1.3": self._SOLEXA13,
            "solexa": self._SOLEXA,
            "sanger": self._SANGER,
        }.get(self.qual_type, self._SANGER)
        return [mapping.get(symbol, 0.0) for symbol in self.quality_string]

    def return_quality_by_pos(self, pos: int) -> float:
        """Retrieve the pre-computed quality score for a single read position."""
        return self.quality_array[pos]


@dataclass
class Protein:
    """Stores the coordinates and description of a protein annotated in the reference."""
    name: str
    description: str
    start_coordinate: int
    end_coordinate: int

    @classmethod
    def from_string(cls, text: str) -> "Protein":
        """Parse a single FASTA header fragment describing a protein interval."""
        left, interval = text.split(":")
        name, desc = left.split("(")
        desc = desc.strip(")")
        start, end = interval.split("-")
        return cls(name=name, description=desc, start_coordinate=int(start), end_coordinate=int(end))


@dataclass
class Gene:
    """Associates a gene name with its protein interval and genomic span."""
    name: str
    start_coordinate: int
    end_coordinate: int
    protein: Protein

    @classmethod
    def from_annotation(cls, line: str) -> "Gene":
        """Build a Gene object from a TSV-style annotation line."""
        parts = line.strip().split("\t")
        name = parts[1]
        start_coordinate = int(parts[2])
        end_coordinate = int(parts[3])
        protein_text = ":".join(parts[4:8])
        protein = Protein.from_string(protein_text)
        return cls(name=name, start_coordinate=start_coordinate, end_coordinate=end_coordinate, protein=protein)


@dataclass
class Amplicon:
    """Captures metadata for each configured amplicon used in the pipeline."""
    label: str
    protein: str
    reference: str
    _5prime_sequence: str
    _3prime_sequence: str
    start_coordinate: int
    end_coordinate: int

    @classmethod
    def from_string(cls, line: str) -> "Amplicon":
        """Translate a CSV/TSV record into the in-memory amplicon definition."""
        clean = line.strip().replace('"', "")
        parts = re.split(r"[\t,\s]+", clean)
        if len(parts) < 7:
            raise ValueError("Amplicon definition requires at least 7 columns")
        return cls(
            label=parts[0],
            protein=parts[1],
            reference=parts[2],
            _5prime_sequence=parts[3],
            _3prime_sequence=parts[4],
            start_coordinate=int(parts[5]),
            end_coordinate=int(parts[6]),
        )


class FullReference:
    """Loads the reference sequence and its annotated proteins from FASTA."""
    def __init__(self, fasta_path: Path) -> None:
        """Initialize the reference by parsing a FASTA file with protein metadata."""
        self.path = fasta_path
        self.seq: str = ""
        self.id: str = ""
        self.proteins: Dict[str, Protein] = {}
        self._load_reference()

    def _load_reference(self) -> None:
        """Populate the sequence and protein dictionary from the FASTA header."""
        record = next(SeqIO.parse(str(self.path), "fasta"))
        self.seq = str(record.seq).upper()
        self.id = record.id
        header_parts = record.description.split(";")
        for part in header_parts:
            part = part.strip()
            if not part:
                continue
            protein = Protein.from_string(part)
            self.proteins[protein.name] = protein

    def get_seq_at(self, position: int, length: int = 3) -> str:
        """Return a substring of the reference centered at the provided base coordinate."""
        return self.seq[position - 1 : position - 1 + length]


GENETIC_CODE = {
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "---": "del",
}

DEFAULT_RATIO_UPPER = 3.16227766
DEFAULT_RATIO_LOWER = 0.316227766
DEFAULT_ENTROPY_THRESHOLD = 0.0


def codon_to_aminoacid(codon: str) -> str:
    """Translate a three-base codon into its corresponding amino acid symbol."""
    return GENETIC_CODE.get(codon.upper(), "X")


@dataclass
class Variant:
    """Tracks counts and strand-aware details for every codon observed at a position."""
    codon: str
    count: int = 0
    fw_reads: int = 0
    rv_reads: int = 0
    balanced_fw: int = 0
    balanced_rv: int = 0
    fw_freq: float = 0.0
    rv_freq: float = 0.0
    ratio: float = 0.0
    amplicons: Dict[str, Dict[str, int]] = field(default_factory=dict)


class SamEntry:
    """Wraps a SAM record with convenience helpers used during variant aggregation."""

    def __init__(self, line: str) -> None:
        """Parse a single SAM line and cache derived attributes that the pipeline needs."""
        self.raw = line.strip()
        self.fields = self.raw.split("\t")
        self.identifier = self.fields[0]
        self.flag = int(self.fields[1])
        self.reference = self.fields[2]
        self.coordinate = int(self.fields[3])
        self.map_quality = int(self.fields[4])
        self.cigar = self.fields[5]
        self.sequence = self.fields[9]
        self.quality_string = self.fields[10]
        self.options = self.fields[11:]
        self.orientation = self._determine_orientation()
        self.occurences = self._parse_occurences()
        self._reference_covered_length = self._compute_reference_covered_length()
        self._ref_to_read = self._build_ref_to_read_map()
        self._quality = Qual(self.quality_string, "sanger")
        self.amplicon = self._parse_amplicon()

    def _determine_orientation(self) -> str:
        """Return whether the read is mapped to the forward or reverse strand."""
        if self.flag & 16:
            return "R"
        return "F"

    def _parse_occurences(self) -> int:
        """Detect if the identifier encodes a multiplicity suffix and return count."""
        parts = self.identifier.split("_")
        tail = parts[-1]
        if tail.isdigit():
            return max(1, int(tail))
        return 1

    def _compute_reference_covered_length(self) -> int:
        """Sum the CIGAR operations that consume reference bases to get coverage."""
        ops = re.findall(r"(\d+)([MIDNSHP=X])", self.cigar)
        total = 0
        for length, op in ops:
            if op in "MDN=X":
                total += int(length)
        return total

    def _build_ref_to_read_map(self) -> Dict[int, Optional[int]]:
        """Walk the CIGAR string and build a mapping from reference positions to read positions.

        CIGAR operations:
        - M/=/X: consume both reference and query
        - I: consume query only (skip in reference mapping)
        - D/N: consume reference only (map to None = gap)
        - S: consume query only (skip in reference mapping)
        - H/P: consume neither
        """
        ref_to_read: Dict[int, Optional[int]] = {}
        ops = re.findall(r"(\d+)([MIDNSHP=X])", self.cigar)
        ref_pos = self.coordinate
        read_pos = 0
        for length_str, op in ops:
            length = int(length_str)
            if op in "M=X":
                for i in range(length):
                    ref_to_read[ref_pos + i] = read_pos + i
                ref_pos += length
                read_pos += length
            elif op == "I":
                read_pos += length
            elif op in "DN":
                for i in range(length):
                    ref_to_read[ref_pos + i] = None
                ref_pos += length
            elif op == "S":
                read_pos += length
            # H and P consume neither
        return ref_to_read

    def _parse_amplicon(self) -> str:
        """Extract the amplicon tag embedded in the read identifier, if present."""
        match = re.search(r"(Amp_[0-9]+)", self.identifier)
        if match:
            return match.group(1)
        return "Amp_NONE"

    def covers(self, position: int) -> bool:
        """Check whether this read spans the requested reference position."""
        return self.coordinate <= position <= self.coordinate + self._reference_covered_length - 1

    def codon_at(self, position: int) -> Optional[str]:
        """Return the codon at the given reference position using the CIGAR-aware mapping."""
        bases = []
        for offset in range(3):
            ref_pos = position + offset
            read_idx = self._ref_to_read.get(ref_pos)
            if read_idx is None and ref_pos in self._ref_to_read:
                bases.append("-")
            elif read_idx is None:
                return None
            else:
                if read_idx < 0 or read_idx >= len(self.sequence):
                    return None
                bases.append(self.sequence[read_idx])
        return "".join(bases)

    def is_mapped(self) -> bool:
        """True when the read is not flagged as unmapped."""
        return not (self.flag & 4)


class SamContainer:
    """Aggregates SAM entries and computes RT variant statistics for an amplicon set."""

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
        self.variants: Dict[int, Dict[str, Variant]] = {}
        self.ratio_upper = ratio_upper
        self.ratio_lower = ratio_lower
        self.entropy_threshold = entropy_threshold

    def load(self) -> None:
        """Load every mapped read from the SAM file so they can be aggregated."""
        with open(self.sam_path, "r") as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                entry = SamEntry(line)
                if entry.is_mapped():
                    self.reads.append(entry)

    def calculate_variant_frequencies(self, protein_name: str, reference: FullReference) -> None:
        """Walk every RT position of the requested protein and build variant histories."""
        if protein_name not in reference.proteins:
            raise ValueError(f"Protein {protein_name} not found in reference")
        protein = reference.proteins[protein_name]
        for pos in range(protein.start_coordinate, protein.end_coordinate, 3):
            variants: Dict[str, Variant] = {}
            fw_depth = 0
            rv_depth = 0
            amplicons_info: Dict[str, Dict[str, int]] = {}
            for read in self.reads:
                # only consider reads that fully span the queried codon
                if not read.covers(pos):
                    continue
                codon = read.codon_at(pos)
                if not codon or len(codon) < 3:
                    continue
                depth_increment = read.occurences
                # track forward/reverse depth separately for ratio calculations
                if read.orientation == "F":
                    fw_depth += depth_increment
                else:
                    rv_depth += depth_increment
                variant = variants.setdefault(codon, Variant(codon=codon))
                variant.count += depth_increment
                if read.orientation == "F":
                    variant.fw_reads += depth_increment
                else:
                    variant.rv_reads += depth_increment
                amplicon_label = read.amplicon or "Amp_NONE"
                # record how each amplicon contributes to the variant counts
                amp_variant = variant.amplicons.setdefault(amplicon_label, {"fw": 0, "rv": 0})
                if read.orientation == "F":
                    amp_variant["fw"] += depth_increment
                else:
                    amp_variant["rv"] += depth_increment
                amp_position = amplicons_info.setdefault(amplicon_label, {"fw_depth": 0, "rv_depth": 0, "depth": 0})
                # accumulate per-amplicon depth metrics for downstream reporting
                if read.orientation == "F":
                    amp_position["fw_depth"] += depth_increment
                else:
                    amp_position["rv_depth"] += depth_increment
                amp_position["depth"] += depth_increment
            entry = {
                "variants": variants,
                "fw_depth": fw_depth,
                "rv_depth": rv_depth,
                "depth": fw_depth + rv_depth,
                "amplicons": amplicons_info,
            }
            # cache the computed metrics for this protein position
            self.variants[pos] = entry
            self.calculate_ratios_nucleotide(pos)
            self.calculate_ratios_aminoacid(pos)
            self.calculate_shannon_entropy_single_position(pos)
            self.log_position_stats(pos)

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

    def write_csv(self, sample_name: str, reference: FullReference, output_path: Path) -> None:
        """Emit the variant summary as a TSV file that mirrors the Perl outputs."""
        fieldnames = [
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
        with open(output_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for pos, entry in sorted(self.variants.items()):
                depth = entry["depth"]
                ratio = entry["fw_depth"] / entry["rv_depth"] if entry["rv_depth"] else 0
                for variant in entry["variants"].values():
                    # emit each variant using the cached depths and ratios
                    writer.writerow(
                        {
                            "FILE": sample_name,
                            "REFERENCE": reference.id,
                            "PROTEIN": "RT",
                            "VARIANT": variant.codon,
                            "POSITION": pos,
                            "FREQ": round((variant.count / depth * 100) if depth else 0, 3),
                            "FWCOV": entry["fw_depth"],
                            "RVCOV": entry["rv_depth"],
                            "TOTALCOV": depth,
                            "RATIO": round(ratio, 3),
                        }
                    )

    def write_xml(self, sample_name: str, reference: FullReference, output_path: Path) -> None:
        """Dump the internal state into XML for downstream diagnostics."""
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        for pos, entry in sorted(self.variants.items()):
            position_elem = ET.SubElement(root, "Position", index=str(pos))
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
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)


def validate_sam_file(path: Path) -> None:
    """Perform a lightweight sanity check on the SAM file before parsing."""
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11 or not fields[0] or not fields[3].isdigit():
                raise ValueError(f"SAM file {path} looks malformed: {line.strip()}")
            return
        raise ValueError(f"SAM file {path} contains no alignment entries")


def validate_reference_file(path: Path) -> None:
    """Ensure the reference FASTA carries the protein metadata that RT parsing needs."""
    record = next(SeqIO.parse(str(path), "fasta"), None)
    if record is None:
        raise ValueError(f"Reference file {path} is empty or not FASTA")
    if "RT" not in record.description:
        raise ValueError(f"Reference {path} header is missing an RT protein annotation")


def validate_amplicon_file(path: Path) -> None:
    """Check that the amplicon TSV/CSV has at least a header plus one valid row."""
    with open(path, "r") as fh:
        header = fh.readline()
        if not header or len(header.strip().split("\t")) < 7:
            raise ValueError(f"Amplicon file {path} does not contain the expected header")
        for line in fh:
            if not line.strip():
                continue
            parts = re.split(r"[\t,\s]+", line.strip().replace('"', ""))
            if len(parts) < 7:
                raise ValueError(f"Amplicon line needs 7 columns: {line.strip()}")
            return
        raise ValueError(f"Amplicon file {path} contains no definitions")


def parse_amplicons(path: Path) -> Dict[str, Amplicon]:
    """Read the amplicon configuration file and construct Amplicon objects."""
    amplicons: Dict[str, Amplicon] = {}
    with open(path, "r") as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            # parse each non-empty line into an Amplicon object
            amplicon = Amplicon.from_string(line)
            amplicons[amplicon.label] = amplicon
    return amplicons


def main() -> None:
    """Command-line entry point that wires inputs/key outputs together."""
    parser = argparse.ArgumentParser(description="Standalone RT variant collector")
    parser.add_argument("sam_file", type=Path)
    parser.add_argument("reference_file", type=Path)
    parser.add_argument("amplicons_file", type=Path)
    parser.add_argument(
        "--ratio-upper",
        type=float,
        default=DEFAULT_RATIO_UPPER,
        help="Upper bound for strand-ratio balancing (default mirrors Perl).",
    )
    parser.add_argument(
        "--ratio-lower",
        type=float,
        default=DEFAULT_RATIO_LOWER,
        help="Lower bound for strand-ratio balancing (default mirrors Perl).",
    )
    parser.add_argument(
        "--entropy-threshold",
        type=float,
        default=DEFAULT_ENTROPY_THRESHOLD,
        help="Minimum Shannon entropy required to mark a position as diverse.",
    )
    args = parser.parse_args()

    if not args.sam_file.exists():
        parser.error(f"SAM file {args.sam_file} does not exist")
    if not args.reference_file.exists():
        parser.error(f"Reference file {args.reference_file} does not exist")
    if not args.amplicons_file.exists():
        parser.error(f"Amplicon file {args.amplicons_file} does not exist")

    logger.info("Validating input formats before parsing")
    validate_sam_file(args.sam_file)
    validate_reference_file(args.reference_file)
    validate_amplicon_file(args.amplicons_file)

    logger.info("Reading amplicon configuration")
    parse_amplicons(args.amplicons_file)

    logger.info("Loading reference")
    # reference contains the RT protein coordinates used for variant calling
    reference = FullReference(args.reference_file)

    logger.info("Parsing SAM entries")
    # load every mapped read from the SAM file into the working buffer
    container = SamContainer(
        args.sam_file,
        ratio_upper=args.ratio_upper,
        ratio_lower=args.ratio_lower,
        entropy_threshold=args.entropy_threshold,
    )
    container.load()

    logger.info("Calculating RT variant frequencies")
    container.calculate_variant_frequencies("RT", reference)

    csv_path = args.sam_file.with_suffix(args.sam_file.suffix + ".tsv")
    xml_path = args.sam_file.with_suffix(args.sam_file.suffix + ".xml")
    logger.info("Writing TSV results")
    # final exports mirror the Perl output format consumed by downstream workflows
    container.write_csv(str(args.sam_file), reference, csv_path)
    logger.info("Writing XML dump")
    container.write_xml(str(args.sam_file), reference, xml_path)


if __name__ == "__main__":
    main()
