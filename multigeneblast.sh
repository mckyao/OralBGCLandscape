#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=48:00:00

MGB_DIR=""          # MultiGeneBlast install directory
GBK_DIR=""          # antiSMASH .gbk files
DB_PATH=""          # database built with makedb
WORK_DIR=""         # per-query MultiGeneBlast output
RESULTS_DIR=""      # collected .txt results
NO_HIT_DIR=""       # .gbk files with no MultiGeneBlast output

MIN_SEQ_COVERAGE=80
MIN_PERCENT_IDENTITY=90
OUT_PAGES=1

export PATH="$MGB_DIR:$PATH"
module load bio/blast-plus

mkdir -p "$WORK_DIR" "$RESULTS_DIR" "$NO_HIT_DIR"

for query in "$GBK_DIR"/*.gbk; do
  [ -e "$query" ] || { echo "No .gbk files in $GBK_DIR"; exit 1; }

  name=$(basename "$query" .gbk)
  out="$WORK_DIR/$(echo "$name" | tr -cd '[:alnum:]')"
  mkdir -p "$out"
  cd "$out"
  export BLASTDB="$out:${BLASTDB:-}"

  len=$(grep "^LOCUS" "$query" | awk '{print $3}')
  if [ -z "$len" ]; then
    echo "Could not read length from LOCUS line: $query"
    continue
  fi

  multigeneblast \
    -in "$query" \
    -db "$DB_PATH" \
    -out "$out" \
    -from 1 \
    -to "$len" \
    -minseqcov "$MIN_SEQ_COVERAGE" \
    -minpercid "$MIN_PERCENT_IDENTITY" \
    -outpages "$OUT_PAGES"

  txt=$(find "$out" -maxdepth 1 -type f -name "*.txt" | head -n 1)
  if [ -n "$txt" ]; then
    mv "$txt" "$RESULTS_DIR/${name}.txt"
  else
    cp "$query" "$NO_HIT_DIR/"
  fi
done
