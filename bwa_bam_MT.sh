#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --array=1-59
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bwamem2_%A_%a.out

set -euo pipefail
module load bio/bwa-mem2
module load bio/samtools

READS_DIR=""       # *.nonhuman.fastq.1.gz / *.nonhuman.fastq.2.gz
OUTPUT_DIR=""
INDEX=""           # indexed BGC reference (.fa or .fna)
BWA_THRESHOLD=80

mkdir -p "${OUTPUT_DIR}/bam" "${OUTPUT_DIR}/stats"

mapfile -t SAMPLES < <(find "${READS_DIR}" -maxdepth 1 -type f -name "*.nonhuman.fastq.1.gz" | sort)
if [ "${SLURM_ARRAY_TASK_ID}" -gt "${#SAMPLES[@]}" ]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds sample count ${#SAMPLES[@]}"
  exit 1
fi

FORWARD="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
REVERSE="${FORWARD/.nonhuman.fastq.1.gz/.nonhuman.fastq.2.gz}"
NAME=$(basename "${FORWARD}" ".nonhuman.fastq.1.gz")
OUT_BAM="${OUTPUT_DIR}/bam/${NAME}.sorted.bam"

[ -f "${REVERSE}" ] || { echo "Missing mate: ${REVERSE}"; exit 1; }
[ -f "${INDEX}" ] || { echo "Missing reference: ${INDEX}"; exit 1; }

bwa-mem2 mem -t "${SLURM_CPUS_PER_TASK}" -T "${BWA_THRESHOLD}" "${INDEX}" "${FORWARD}" "${REVERSE}" \
  | samtools view -@ "${SLURM_CPUS_PER_TASK}" -b - \
  | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${OUT_BAM}" -

samtools index -@ "${SLURM_CPUS_PER_TASK}" "${OUT_BAM}"
samtools flagstat -@ "${SLURM_CPUS_PER_TASK}" "${OUT_BAM}" > "${OUTPUT_DIR}/stats/${NAME}.flagstat.txt"
