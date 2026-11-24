#!/usr/bin/env python3
"""Count reads per reference from BWA BAM/SAM files.

Reads that map equally well to several references (via XA:Z) are split
fractionally (1/N) across them.
"""
import argparse
import os
import re
import subprocess
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", required=True, nargs="+", help="BAM/SAM files")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def alt_hits(edit, cols):
    hits = []
    for field in cols[11:]:
        if field.startswith("XA:Z:"):
            for alt in field[5:].split(";"):
                if not alt:
                    continue
                parts = alt.split(",")
                if parts[-1].isdigit() and int(parts[-1]) <= edit:
                    hits.append(parts[0])
    return hits


def parse_line(line, hits_by_read):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 11:
        return
    edit = None
    for field in cols[11:]:
        if field.startswith("NM:i:"):
            try:
                edit = int(field.split(":")[2])
            except ValueError:
                edit = None
    if edit is None:
        return
    refs = hits_by_read.setdefault(cols[0], [])
    refs.append(cols[2])
    refs.extend(alt_hits(edit, cols))


def output_name(path):
    base = os.path.basename(path)
    new = re.sub(r"\.(sam|bam)$", ".txt", base, flags=re.IGNORECASE)
    return new if new != base else base + ".txt"


def main():
    args = get_args()
    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, output_name(args.files[0]))

    counts = {}
    hits_by_read = {}

    for path in args.files:
        if args.verbose:
            print(f"Loading {path}", file=sys.stderr)
        # Drop unmapped (0x4) and supplementary (0x800) alignments
        cmd = f"samtools view -h -F 0x4 {path} | samtools view -h -F 0x800 -"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        for line in proc.stdout:
            if line.startswith("@SQ"):
                for field in line.rstrip("\n").split("\t"):
                    if field.startswith("SN:"):
                        counts.setdefault(field[3:], 0.0)
                continue
            if line.startswith("@"):
                continue
            parse_line(line, hits_by_read)
        proc.stdout.close()
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"samtools failed for {path} (exit {proc.returncode})")

    for refs in hits_by_read.values():
        unique = set(refs)
        if not unique:
            continue
        inc = 1.0 / len(unique)
        for r in unique:
            counts[r] = counts.get(r, 0.0) + inc

    with open(output_file, "w") as f:
        for r in sorted(counts):
            f.write(f"{r}\t{counts[r]}\n")

    if args.verbose:
        print(f"Wrote {output_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
