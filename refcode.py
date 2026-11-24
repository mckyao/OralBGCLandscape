#!/usr/bin/env python3
"""Merge per-sample BGC count files into one abundance table."""
import argparse
import glob
import os
import re
import pandas as pd


def sample_name(path):
    name = os.path.basename(path)
    name = re.sub(r"\.txt$", "", name)
    name = re.sub(r"_result$", "", name)
    name = re.sub(r"\.sorted$", "", name)
    return name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_folder", required=True, help="Folder of per-sample count files")
    parser.add_argument("-output_path", default="bgc_abundance_table_final.txt")
    parser.add_argument("-pattern", default="*result.txt", help="Glob for count files")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_folder, args.pattern)))
    if not files:
        raise SystemExit(f"No files matching {args.pattern} in {args.input_folder}")

    line_counts = {sum(1 for _ in open(f)) for f in files}
    if len(line_counts) != 1:
        print("Warning: count files have different numbers of lines.")

    with open(files[0]) as f:
        rownames = [line.split()[0] for line in f]
    df = pd.DataFrame(index=rownames)

    for path in files:
        with open(path) as f:
            df[sample_name(path)] = [line.split()[1] for line in f]

    df = df[sorted(df.columns)]
    df.to_csv(args.output_path, sep="\t", header=True, index=True)
    print(f"Wrote {args.output_path}")


if __name__ == "__main__":
    main()
