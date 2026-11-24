#!/bin/bash
#SBATCH --job-name=
#SBATCH --account=
#SBATCH --partition=
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-40
#SBATCH --output=logs/trim_%A_%a.out

set -euo pipefail
module load java

RAW_DIR=""          # raw *_1.fastq / *_2.fastq
TRIM_DIR=""         # output directory
TRIMMOMATIC_JAR=""  # path to trimmomatic-0.39.jar
ADAPTERS=""         # path to TruSeq3-PE.fa
SRR_LIST="srr_list.txt"  # one accession per line

mkdir -p "${TRIM_DIR}"
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SRR_LIST}")

java -jar "${TRIMMOMATIC_JAR}" PE \
  -threads "${SLURM_CPUS_PER_TASK}" \
  "${RAW_DIR}/${ACCESSION}_1.fastq" "${RAW_DIR}/${ACCESSION}_2.fastq" \
  "${TRIM_DIR}/${ACCESSION}_1_paired_trimmed.fastq" "${TRIM_DIR}/${ACCESSION}_1_unpaired_trimmed.fastq" \
  "${TRIM_DIR}/${ACCESSION}_2_paired_trimmed.fastq" "${TRIM_DIR}/${ACCESSION}_2_unpaired_trimmed.fastq" \
  ILLUMINACLIP:"${ADAPTERS}":2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 MINLEN:36

echo "Finished ${ACCESSION}"
