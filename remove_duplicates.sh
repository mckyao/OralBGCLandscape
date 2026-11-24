#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=05:00:00

module load python

MGB_RESULTS=""   # MultiGeneBlast .txt results
GBK_DIR=""       # antiSMASH .gbk files
OUT_DIR=""       # non-redundant .gbk output

python3 remove_duplicates.py --mgb "$MGB_RESULTS" --gbk "$GBK_DIR" --out "$OUT_DIR"
