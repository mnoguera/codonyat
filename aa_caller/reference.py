from __future__ import annotations

from pathlib import Path
from collections.abc import Sequence
from typing import Dict, List

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
        record = next(SeqIO.parse(str(self.path), "fasta"), None)
        if record is None:
            raise ValueError(f"Reference file {self.path} is empty or not valid FASTA")
        self.seq = str(record.seq).upper()
        self.id = record.id
        header_parts = record.description.split(";")
        for i, part in enumerate(header_parts):
            part = part.strip()
            if not part:
                continue
            # The first segment includes the sequence ID prefix (e.g.
            # "K03455|HIVHXB2CG RT(Protease):1-99"); strip it so
            # Protein.from_string receives just "RT(Protease):1-99".
            if i == 0 and part.startswith(record.id):
                part = part[len(record.id):].strip()
                if not part:
                    continue
            protein = Protein.from_string(part)
            self.proteins[protein.name] = protein

    def get_seq_at(self, position: int, length: int = 3) -> str:
        """Return a substring of the reference centered at the provided base coordinate."""
        return self.seq[position - 1 : position - 1 + length]


def parse_protein_spec(spec: str | Sequence[str]) -> List[str]:
    """Split a protein selection (``"RT"``, ``"PR,RT,INT"``, ``"all"`` or a list) into names.

    Order is kept and duplicates are dropped. ``"all"`` is returned as ``["all"]``
    and expanded by :func:`resolve_proteins` once the reference is loaded.
    """
    if isinstance(spec, str):
        items = [part.strip() for part in spec.split(",")]
    else:
        items = [str(part).strip() for part in spec]
    names: List[str] = []
    for item in items:
        if item and item not in names:
            names.append(item)
    if not names:
        raise ValueError("No protein selected")
    if "all" in (name.lower() for name in names) and len(names) > 1:
        raise ValueError("'all' cannot be combined with other protein names")
    return names


def resolve_proteins(spec: str | Sequence[str], reference: FullReference) -> List[str]:
    """Resolve a protein selection against the reference annotations.

    ``"all"`` selects every annotated protein in FASTA-header order.
    """
    names = parse_protein_spec(spec)
    if len(names) == 1 and names[0].lower() == "all":
        if not reference.proteins:
            raise ValueError(f"Reference {reference.path} has no annotated proteins")
        return list(reference.proteins)
    missing = [name for name in names if name not in reference.proteins]
    if missing:
        raise ValueError(f"Protein(s) not found in reference: {', '.join(missing)}")
    return names
