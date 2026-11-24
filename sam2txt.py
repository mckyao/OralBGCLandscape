#!/usr/bin/env python3
"""Count reads per reference from BWA SAM/BAM files, with alignment filters."""
import argparse
import os
import re
import subprocess


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", required=True, nargs="+", help="SAM/BAM files")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-r", "--rate", action="store_true",
                        help="Split multi-mapped reads fractionally (1/N)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--min_len", type=int, default=1500, help="Min alignment length (bp)")
    parser.add_argument("--min_identity", type=float, default=0.9, help="Min identity (0-1)")
    parser.add_argument("--min_coverage", type=float, default=0.7,
                        help="Min fraction of reference covered")
    return parser.parse_args()


def aligned_length(cigar):
    return sum(int(n) for n, _ in re.findall(r"(\d+)([M=X])", cigar))


def nm_and_xa(cols):
    nm, xa = None, None
    for c in cols[11:]:
        if c.startswith("NM:i:"):
            try:
                nm = int(c.split(":")[-1])
            except ValueError:
                nm = None
        elif c.startswith("XA:Z:"):
            xa = c[5:]
    return nm, xa


def passes(aln_len, mismatches, ref_len, min_len, min_id, min_cov):
    if aln_len < min_len or aln_len == 0 or not ref_len:
        return False
    return (aln_len - mismatches) / aln_len >= min_id and aln_len / ref_len >= min_cov


def alt_hits(xa, primary_nm, ref_lengths, min_len, min_id, min_cov):
    if xa is None:
        return []
    hits = []
    # XA format: ref,pos,strand,CIGAR,NM;...
    for entry in xa.split(";"):
        fields = entry.split(",")
        if len(fields) < 5:
            continue
        try:
            mm = int(fields[4])
        except ValueError:
            continue
        if primary_nm is not None and mm > primary_nm:
            continue
        ref = fields[0]
        if passes(aligned_length(fields[3]), mm, ref_lengths.get(ref), min_len, min_id, min_cov):
            hits.append(ref)
    return hits


def parse_line(line, hits_by_read, ref_lengths, min_len, min_id, min_cov):
    cols = line.split()
    if len(cols) < 11:
        return
    ref = cols[2]
    aln_len = aligned_length(cols[5])
    nm, xa = nm_and_xa(cols)
    if nm is None:
        return
    if not passes(aln_len, nm, ref_lengths.get(ref), min_len, min_id, min_cov):
        return
    refs = hits_by_read.setdefault(cols[0], [])
    refs.append(ref)
    refs.extend(alt_hits(xa, nm, ref_lengths, min_len, min_id, min_cov))


def main():
    args = get_args()
    os.makedirs(args.output, exist_ok=True)

    base = os.path.basename(args.files[0])
    out_name = re.sub(r"\.(sam|bam)$", ".txt", base)
    if out_name == base:
        out_name = base + ".txt"
    output_file = os.path.join(args.output, out_name)

    counts, ref_lengths, hits_by_read = {}, {}, {}

    for path in args.files:
        proc = subprocess.Popen(f"samtools view -H {path}", shell=True,
                                stdout=subprocess.PIPE, universal_newlines=True)
        for line in proc.stdout:
            if not line.startswith("@SQ"):
                continue
            sn = ln = None
            for p in line.split():
                if p.startswith("SN:"):
                    sn = p[3:]
                elif p.startswith("LN:"):
                    ln = p[3:]
            if sn and ln:
                ref_lengths[sn] = int(ln)
                counts.setdefault(sn, 0)

    for path in args.files:
        if args.verbose:
            print(f"Processing {path}")
        # Drop unmapped (0x4) and supplementary (0x800) alignments
        cmd = f"samtools view -h -F 0x4 {path} | samtools view -h -F 0x800 -"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        for line in proc.stdout:
            if not line.startswith("@"):
                parse_line(line, hits_by_read, ref_lengths,
                           args.min_len, args.min_identity, args.min_coverage)

    for refs in hits_by_read.values():
        unique = set(refs)
        if not unique:
            continue
        inc = 1.0 / len(unique) if args.rate else 1.0
        for r in unique:
            counts[r] = counts.get(r, 0) + inc

    with open(output_file, "w") as f:
        for r in sorted(counts):
            f.write(f"{r}\t{counts[r]}\n")

    if args.verbose:
        print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
