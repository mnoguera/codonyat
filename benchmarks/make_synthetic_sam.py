#!/usr/bin/env python3
"""Generate a deterministic synthetic SAM file for codonyat benchmarks and regression tests.

Reads are sampled from an annotated reference (e.g. HXB2 pol) and carry
substitutions, small insertions/deletions, soft clips, both strands, a few
unmapped records, ``Amp_<n>`` tags and occasional ``_<n>`` multiplicity
suffixes, so every code path of the codon counter is exercised.

Usage:
    python benchmarks/make_synthetic_sam.py REFERENCE.fasta OUT.sam --reads 1000000 \
        [--start 2253] [--end 5093] [--read-length 250] [--seed 1]
"""

from __future__ import annotations

import argparse
import random


def read_fasta(path: str) -> tuple[str, str]:
    name, seq = "", []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip().upper())
    return name, "".join(seq)


def mutate(rng: random.Random, seg: str) -> tuple[str, str]:
    """Return (read_sequence, cigar_core) for a reference segment."""
    bases = "ACGT"
    out: list[str] = []
    ops: list[tuple[str, int]] = []

    def push(op: str, n: int = 1) -> None:
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + n)
        else:
            ops.append((op, n))

    i = 0
    while i < len(seg):
        r = rng.random()
        if r < 0.0015 and 5 < i < len(seg) - 5:  # deletion of 1 or 3 bases
            n = rng.choice((1, 3))
            push("D", n)
            i += n
            continue
        if r < 0.003 and 5 < i < len(seg) - 5:  # insertion of 1 or 3 bases
            n = rng.choice((1, 3))
            out.extend(rng.choice(bases) for _ in range(n))
            push("I", n)
        b = seg[i]
        if rng.random() < 0.01:
            b = rng.choice(bases.replace(b, "") or bases)
        elif rng.random() < 0.0005:
            b = "N"
        out.append(b)
        push("M")
        i += 1
    return "".join(out), "".join(f"{n}{op}" for op, n in ops)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("reference")
    ap.add_argument("out")
    ap.add_argument("--reads", type=int, default=100_000)
    ap.add_argument("--start", type=int, default=2253)
    ap.add_argument("--end", type=int, default=5093)
    ap.add_argument("--read-length", type=int, default=250)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref_name, ref = read_fasta(args.reference)
    # A few common "haplotype" codon changes so positions carry real mixtures.
    hot = {p: rng.choice("ACGT") for p in range(args.start, args.end, 37)}
    qual_chars = "#+5:?FI"
    with open(args.out, "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{ref_name}\tLN:{len(ref)}\n")
        for n in range(args.reads):
            length = max(30, int(rng.gauss(args.read_length, 25)))
            pos = rng.randint(max(1, args.start - 100), max(args.start, args.end - length // 2))
            seg = list(ref[pos - 1 : pos - 1 + length])
            for p, b in hot.items():
                if pos <= p < pos + len(seg) and rng.random() < 0.3:
                    seg[p - pos] = b
            seq, core = mutate(rng, "".join(seg))
            lead = rng.randint(0, 8) if rng.random() < 0.3 else 0
            trail = rng.randint(0, 8) if rng.random() < 0.3 else 0
            seq = "".join(rng.choice("ACGT") for _ in range(lead)) + seq + "".join(rng.choice("ACGT") for _ in range(trail))
            cigar = (f"{lead}S" if lead else "") + core + (f"{trail}S" if trail else "")
            flag = 16 if rng.random() < 0.5 else 0
            if rng.random() < 0.01:
                flag |= 4
            amp = f"Amp_{rng.randint(1, 3)}" if rng.random() < 0.8 else "noamp"
            qname = f"r{n}:{amp}"
            if rng.random() < 0.05:
                qname += f"_{rng.randint(1, 5)}"
            qual = "".join(rng.choice(qual_chars) for _ in range(len(seq)))
            fh.write(f"{qname}\t{flag}\t{ref_name}\t{pos}\t{rng.randint(20, 42)}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")


if __name__ == "__main__":
    main()
