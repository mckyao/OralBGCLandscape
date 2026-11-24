#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Usage: sbatch run_bgc_blast.sh <bgc.fa> [read_len=150] [pct_id=97]

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch $0 <bgc_fasta> [read_len=150] [pct_id=97]" >&2
  exit 1
fi

BGC_FASTA="$1"
READ_LEN="${2:-150}"
PCT_ID="${3:-97}"
PREFIX="$(basename "${BGC_FASTA%.*}")"

module load bio/blast-plus

makeblastdb -in "$BGC_FASTA" -dbtype nucl -out "${PREFIX}_db"

BLAST_OUT="${PREFIX}_all_vs_all.blast.tsv"
blastn \
  -query "$BGC_FASTA" \
  -db "${PREFIX}_db" \
  -out "$BLAST_OUT" \
  -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send" \
  -num_threads "$SLURM_CPUS_PER_TASK" \
  -max_target_seqs 100000 \
  -perc_identity "$PCT_ID"

# Keep non-self hits at or above the length and identity thresholds
awk -v min_len="$READ_LEN" -v min_pid="$PCT_ID" \
  '$1 != $2 && $3 >= min_pid && $4 >= min_len' \
  "$BLAST_OUT" > "${PREFIX}_ambiguous_windows.len${READ_LEN}.pid${PCT_ID}.tsv"
