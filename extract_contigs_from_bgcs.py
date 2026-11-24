#!/usr/bin/env python3
"""Extract the assembly contigs that each BGC was called from.

Usage:
  python extract_contigs_from_bgcs.py \
      --bgc-fasta bgcs.fa \
      --assembly-dir assemblies/ \
      --output-fasta bgc_contigs.fa
"""
import argparse
import os
import re
import sys

REGION_SUFFIX = re.compile(r"(\.region\d+|_region\d+|\.cluster\d+|_c\d{5})$")
FASTA_EXTS = (".fa", ".fasta", ".fna")


def fasta_iter(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def strip_region_suffix(token):
    m = REGION_SUFFIX.search(token)
    return token[:m.start()] if m else token


def index_assemblies(assembly_dir):
    # header token -> set of assembly files containing it
    index = {}
    for root, _, files in os.walk(assembly_dir):
        for fn in files:
            if not fn.lower().endswith(FASTA_EXTS):
                continue
            path = os.path.join(root, fn)
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        for tok in line[1:].split():
                            index.setdefault(tok, set()).add(path)
    return index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bgc-fasta", required=True, help="antiSMASH BGC multi-FASTA")
    ap.add_argument("--assembly-dir", required=True, help="Directory of assembly FASTA files")
    ap.add_argument("--output-fasta", required=True, help="Output FASTA of matched contigs")
    args = ap.parse_args()

    index = index_assemblies(args.assembly_dir)
    if not index:
        sys.exit(f"No FASTA files found in {args.assembly_dir}")

    file_to_tokens = {}
    n_bgcs = n_matched = 0
    for header, _ in fasta_iter(args.bgc_fasta):
        n_bgcs += 1
        # Try the first header token, then the second
        for tok in [strip_region_suffix(t) for t in header.split()[:2]]:
            if tok in index:
                n_matched += 1
                for path in index[tok]:
                    file_to_tokens.setdefault(path, set()).add(tok)
                break

    if not file_to_tokens:
        sys.exit("No matching contigs found.")

    written = set()
    with open(args.output_fasta, "w") as out:
        for path, tokens in file_to_tokens.items():
            for header, seq in fasta_iter(path):
                hit = next((t for t in header.split() if t in tokens), None)
                if hit and (path, hit) not in written:
                    written.add((path, hit))
                    out.write(f">{header}\n{seq}\n")

    print(f"BGCs: {n_bgcs}, matched: {n_matched}, contigs written: {len(written)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
