from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .constants import DEFAULT_ENTROPY_THRESHOLD, DEFAULT_RATIO_LOWER, DEFAULT_RATIO_UPPER
from .container import SamContainer
from .reference import FullReference
from .validators import parse_amplicons, validate_amplicon_file, validate_reference_file, validate_sam_file

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


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
