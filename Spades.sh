#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=0-80
#SBATCH --output=logs/spades_%A_%a.out

set -euo pipefail
module load anaconda3
source activate spades

INPUT_DIR=""         # trimmed paired reads
WORK_DIR=""          # per-sample SPAdes output
SCAFFOLD_DIR=""      # final renamed scaffolds

mkdir -p "$WORK_DIR" "$SCAFFOLD_DIR"

files=("$INPUT_DIR"/*_1_paired_trimmed.fastq)
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#files[@]}" ]; then
  echo "Array index $SLURM_ARRAY_TASK_ID exceeds file count ${#files[@]}"
  exit 1
fi

forward="${files[$SLURM_ARRAY_TASK_ID]}"
accession=$(basename "$forward" "_1_paired_trimmed.fastq")
reverse="$INPUT_DIR/${accession}_2_paired_trimmed.fastq"

if [[ ! -f "$reverse" ]]; then
  echo "Missing reverse read for $accession"
  exit 1
fi

out="$WORK_DIR/$accession"
mkdir -p "$out"

spades.py --meta -t "$SLURM_CPUS_PER_TASK" -m 350 -1 "$forward" -2 "$reverse" -o "$out"

if [[ -f "$out/scaffolds.fasta" ]]; then
  mv "$out/scaffolds.fasta" "$SCAFFOLD_DIR/${accession}_scaffolds.fasta"
else
  echo "scaffolds.fasta not found for $accession"
fi
