#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G

MMSEQS_BIN=""     # directory containing the mmseqs binary
QUERY_FASTA=""    # contigs from extract_contigs_from_bgcs.py
GTDB_DB=""        # MMseqs2 GTDB database
OUT_PREFIX=""

export PATH="$MMSEQS_BIN:$PATH"

mmseqs easy-taxonomy "$QUERY_FASTA" "$GTDB_DB" "$OUT_PREFIX" tmp \
  --threads "$SLURM_CPUS_PER_TASK" \
  -s 5 \
  --max-seqs 300 \
  --lca-mode 3 \
  --tax-lineage 1 \
  --blacklist "" \
  --orf-filter 0
