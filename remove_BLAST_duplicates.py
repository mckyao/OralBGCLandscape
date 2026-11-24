#!/usr/bin/env python3
"""Delete BGC .fasta/.gbk files whose core ID appears in a removal list."""
import argparse
import os
from collections import defaultdict

VALID_EXTS = {".fasta", ".fa", ".fna", ".gbk", ".gbff"}


def extract_core_id(name):
    # NODE_1_length_2758_cov_3.01.region001.gbk -> NODE_1_length_2758
    base = os.path.splitext(name)[0]
    base = base.split(".region")[0]
    base = base.split("_cov")[0]
    return base.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--list", required=True, help="BGC IDs to remove, one per line")
    parser.add_argument("--folder", required=True, help="Folder of .fasta/.gbk files")
    args = parser.parse_args()

    with open(args.list) as f:
        remove_ids = {extract_core_id(line.strip()) for line in f if line.strip()}

    core_to_files = defaultdict(list)
    for filename in os.listdir(args.folder):
        if os.path.splitext(filename)[1].lower() in VALID_EXTS:
            core_to_files[extract_core_id(filename)].append(filename)

    removed = 0
    missing = 0
    for core_id in sorted(remove_ids):
        if core_id not in core_to_files:
            print(f"No files found for {core_id}")
            missing += 1
            continue
        for filename in core_to_files[core_id]:
            os.remove(os.path.join(args.folder, filename))
            removed += 1

    print(f"Removed {removed} files. {missing} IDs had no matching files.")


if __name__ == "__main__":
    main()
