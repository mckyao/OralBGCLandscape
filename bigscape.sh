#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=24:00:00

module load anaconda3
source activate bigscape

GBK_DIR=""     # non-redundant .gbk files
OUT_DIR=""
PFAM=""        # Pfam-A.hmm

# Add --mibig-version 4.0 to include MIBiG reference BGCs
bigscape cluster \
  -i "$GBK_DIR" \
  -o "$OUT_DIR" \
  -p "$PFAM" \
  -c "$SLURM_CPUS_PER_TASK" \
  --alignment-mode auto \
  --include-singletons \
  --gcf-cutoffs 0.6 \
  --mix
