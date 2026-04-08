from __future__ import annotations

from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .models import Protein


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
