from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List


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
