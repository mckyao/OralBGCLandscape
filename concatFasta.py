#!/usr/bin/env python3
"""Concatenate all .fasta files in a directory into one file."""
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Directory of .fasta files")
    parser.add_argument("--output", required=True, help="Combined output (.fa)")
    args = parser.parse_args()

    with open(args.output, "w") as out:
        for filename in sorted(os.listdir(args.input_dir)):
            if filename.endswith(".fasta"):
                with open(os.path.join(args.input_dir, filename)) as f:
                    out.write(f.read())


if __name__ == "__main__":
    main()
