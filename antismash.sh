#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --array=1-42
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --output=logs/antismash_%A_%a.out

set -euo pipefail
module load python
conda activate antismash

SCAFFOLD_DIR=""   # SPAdes scaffolds (*.fasta)
OUT_BASE=""       # antiSMASH output

mkdir -p logs "$OUT_BASE"

FILES=("$SCAFFOLD_DIR"/*.fasta)
if [[ $SLURM_ARRAY_TASK_ID -gt ${#FILES[@]} ]]; then
  echo "Array index $SLURM_ARRAY_TASK_ID exceeds file count ${#FILES[@]}"
  exit 1
fi

FASTA_FILE=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$(basename "$FASTA_FILE" .fasta)

antismash "$FASTA_FILE" \
  --output-dir "$OUT_BASE/$NAME" \
  --cpu "$SLURM_CPUS_PER_TASK" \
  --cb-knownclusters --tfbs --cb-subclusters --asf --rre \
  --genefinding-tool prodigal
