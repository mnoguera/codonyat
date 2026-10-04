"""Light wrapper to expose the codon-aware variant caller as a package."""

from .cli import main
from .container import SamContainer
from .reference import FullReference, parse_protein_spec, resolve_proteins
from .runner import VariantCallResult, call_variants, call_variants_from_args, runner_cli
from .sam import SamEntry, iter_sam_lines
from .validators import parse_amplicons

__all__ = [
    "call_variants",
    "call_variants_from_args",
    "runner_cli",
    "VariantCallResult",
    "main",
    "FullReference",
    "SamContainer",
    "SamEntry",
    "parse_amplicons",
    "parse_protein_spec",
    "resolve_proteins",
    "iter_sam_lines",
]
