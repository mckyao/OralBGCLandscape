#!/usr/bin/env python3
"""Convert each .gbk file in a directory to .fasta."""
import argparse
import os
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Directory of .gbk files")
    parser.add_argument("--output-dir", required=True, help="Directory for .fasta files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for filename in os.listdir(args.input_dir):
        if not filename.endswith(".gbk"):
            continue
        in_path = os.path.join(args.input_dir, filename)
        out_path = os.path.join(args.output_dir, os.path.splitext(filename)[0] + ".fasta")
        with open(in_path) as fin, open(out_path, "w") as fout:
            for rec in SeqIO.parse(fin, "genbank"):
                fout.write(f">{rec.id} {rec.description}\n{rec.seq}\n")


if __name__ == "__main__":
    main()
