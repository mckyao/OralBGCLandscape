#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --time=24:00:00

set -euo pipefail
module load python
module load bio/samtools/1.17-gcc-11.4.0

SCRIPT_DIR=""                          # directory containing bam2txt.py and refcode.py
BAM_DIR=""                             # sorted BAM files from bwa_bam_MT.sh
TXT_DIR="${BAM_DIR}/txt"
COUNTS_TABLE="${BAM_DIR}/bgc_counts.txt"

mkdir -p "${TXT_DIR}"

for bam in "${BAM_DIR}"/*.bam; do
  python3 "${SCRIPT_DIR}/bam2txt.py" -f "${bam}" -o "${TXT_DIR}"
done

python3 "${SCRIPT_DIR}/refcode.py" \
  -input_folder "${TXT_DIR}" \
  -output_path "${COUNTS_TABLE}" \
  -pattern "*.txt"
