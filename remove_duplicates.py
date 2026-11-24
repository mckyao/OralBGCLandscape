#!/usr/bin/env python3
"""Pick non-redundant BGCs from MultiGeneBlast results.

Pass 1: drop a subject BGC if it has the same gene count as the query and
        every query gene hits it at >90% identity.
Pass 2: drop a query BGC if >=90% of its genes hit a larger subject BGC at
        >90% identity.
Remaining BGCs are copied to the output folder.
"""
import argparse
import glob
import os
import re
import shutil

IDENTITY_CUTOFF = 90.0
HIT_RATIO_CUTOFF = 0.9


def bgc_name(path):
    return os.path.splitext(os.path.basename(path))[0].strip()


def count_query_genes(lines):
    total, reading = 0, False
    for line in lines:
        line = line.strip()
        if re.search(r"Table of genes,\s*locations,\s*strands and annotations of query cluster", line, re.I):
            reading = True
            continue
        if reading:
            if line == "" or re.match(r"Table of", line, re.I):
                reading = False
            else:
                total += 1
    return total


def parse_subject_blocks(lines):
    subjects = []
    subject, gene_count, identities = None, 0, []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()

        m = re.match(r"^\d+\.\s*(\S+)", line)
        if m:
            subject, gene_count, identities = m.group(1), 0, []
            i += 1
            continue

        if re.search(r"Table of genes,\s*locations,\s*strands and annotations of subject cluster", line, re.I):
            i += 1
            while i < n:
                curr = lines[i].strip()
                if curr == "" or re.match(r"Table of", curr, re.I):
                    break
                gene_count += 1
                i += 1
            continue

        if re.search(r"Table of Blast hits", line, re.I):
            i += 1
            while i < n:
                curr = lines[i].strip()
                if curr == "" or re.match(r"Table of", curr, re.I):
                    break
                if not (curr.startswith("Query") or curr.startswith("-")):
                    fields = curr.split()
                    if len(fields) >= 6:
                        try:
                            identities.append(float(fields[2]))
                        except ValueError:
                            identities.append(0.0)
                i += 1
            if subject:
                subjects.append({
                    "subject": subject.strip(),
                    "subject_gene_count": gene_count,
                    "blast_hits_count": len(identities),
                    "identities": identities,
                })
            continue
        i += 1
    return subjects


def all_high(identities):
    return bool(identities) and all(x > IDENTITY_CUTOFF for x in identities)


def load(path):
    with open(path) as f:
        lines = f.readlines()
    return count_query_genes(lines), parse_subject_blocks(lines)


def exact_duplicate(path):
    query = bgc_name(path)
    n_genes, blocks = load(path)
    for b in blocks:
        if (b["blast_hits_count"] == n_genes and b["subject_gene_count"] == n_genes
                and all_high(b["identities"]) and b["subject"] != query):
            return b["subject"]
    return None


def contained_in_larger(path):
    n_genes, blocks = load(path)
    for b in blocks:
        ratio = b["blast_hits_count"] / n_genes if n_genes else 0
        if (ratio >= HIT_RATIO_CUTOFF and all_high(b["identities"])
                and b["subject_gene_count"] > b["blast_hits_count"]):
            return n_genes, dict(b, ratio=ratio)
    return n_genes, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mgb", required=True, help="MultiGeneBlast .txt results")
    ap.add_argument("--gbk", required=True, help="Original .gbk files")
    ap.add_argument("--out", required=True, help="Output folder for non-redundant .gbk")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    mgb_files = glob.glob(os.path.join(args.mgb, "*.txt"))
    file_by_name = {bgc_name(f): f for f in mgb_files}
    candidates = list(file_by_name)

    exact = []
    for name in list(candidates):
        dup = exact_duplicate(file_by_name[name])
        if dup:
            exact.append((name, dup, file_by_name[name]))
            if dup in candidates:
                candidates.remove(dup)

    with open(os.path.join(args.mgb, "exact_duplicates_exclusions.txt"), "w") as f:
        for q, d, src in exact:
            f.write(f"Query: {q}\tDuplicate: {d}\tFile: {src}\n")

    summary = []
    for name in candidates:
        n_genes, hit = contained_in_larger(file_by_name[name])
        summary.append((name, n_genes, hit))

    with open(os.path.join(args.out, "master_parent_groups.txt"), "w") as f:
        f.write("=== Exact duplicates ===\n")
        for q, d, src in exact:
            f.write(f"Query: {q}\tDuplicate: {d}\tFile: {src}\n")
        if not exact:
            f.write("None\n")
        f.write("\n=== Candidate evaluations ===\n")
        for i, (name, n_genes, hit) in enumerate(summary, 1):
            f.write(f"Group {i}: {name}\nGenes: {n_genes}\n")
            if hit:
                f.write("EXCLUDED\n")
                f.write(f"  Subject: {hit['subject']}\n")
                f.write(f"  Blast hits: {hit['blast_hits_count']} (ratio {hit['ratio']:.2f})\n")
                f.write(f"  Subject genes: {hit['subject_gene_count']}\n")
                f.write(f"  Identities: {hit['identities']}\n")
            else:
                f.write("INCLUDED\n")
            f.write("-" * 40 + "\n")

    for name, _, hit in summary:
        if hit:
            continue
        matches = glob.glob(os.path.join(args.gbk, name + "*.gbk"))
        if matches:
            shutil.copy(matches[0], os.path.join(args.out, name + ".gbk"))
        else:
            print(f"No .gbk found for {name}")


if __name__ == "__main__":
    main()
