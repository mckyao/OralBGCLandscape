# BGCLandscape Pipeline

Pipeline for detecting, deduplicating, quantifying and comparing biosynthetic gene clusters (BGCs) in oral metagenomes. The scripts are written for SLURM; set `--job-name`, `--account`, `--partition` and the module names for your system.

Each script has empty path variables at the top. Fill these in before running.

## Workflow

```
SRA download → Trimmomatic → metaSPAdes → antiSMASH
  → MultiGeneBlast dedup → BLAST dedup → BWA-MEM2 mapping
  → BiG-SCAPE → R statistics → Cytoscape → MMseqs2 taxonomy
```

## Requirements

SRA Tools, Trimmomatic 0.39, SPAdes, antiSMASH, MultiGeneBlast, BLAST+, BWA-MEM2, samtools, BiG-SCAPE, MMseqs2 with a GTDB database, Pfam-A.hmm, Cytoscape.

Python 3: `biopython`, `pandas`
R: `dplyr`, `stringr`

## Steps

### 1. Download reads

```bash
fasterq-dump <SRR_ID> --outdir raw_reads/ --threads 8 --split-files
```

Put one accession per line in `srr_list.txt`.

### 2. Trim reads

```bash
sbatch Trimmomatic.sh
```

Settings: `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36`. Only the paired outputs are used downstream.

### 3. Assemble

```bash
sbatch Spades.sh
```

Runs `spades.py --meta -m 350` per sample. Outputs `<accession>_scaffolds.fasta`.

### 4. Detect BGCs

```bash
sbatch antismash.sh
```

Then collect all region `.gbk` files into one directory:

```bash
find <antismash_output> -type f -name "*.region*.gbk" -exec mv {} <all_gbk_dir>/ \;
```

### 5. Deduplicate, round 1 (MultiGeneBlast)

```bash
makedb <db_name> <all_gbk_dir>
sbatch multigeneblast.sh
sbatch remove_duplicates.sh
```

`remove_duplicates.py` drops BGCs that are exact duplicates or where ≥90% of genes hit a larger BGC at >90% identity.

### 6. Deduplicate, round 2 (BLAST)

Convert the round-1 `.gbk` set to one FASTA:

```bash
python3 gbk2fasta.py --input-dir <gbk_dir> --output-dir <fasta_dir>
python3 concatFasta.py --input-dir <fasta_dir> --output bgcs.fa
```

Run all-vs-all BLAST (args: FASTA, read length, percent identity):

```bash
sbatch run_bgc_blast.sh bgcs.fa 150 96
```

Set `input_file` in `BGCcoverage.r` to the `*_ambiguous_windows.*.tsv` output and run it. It writes `bgcs_to_remove.txt`, which lists BGCs that are ≥95% covered by a longer BGC.

```bash
python3 remove_BLAST_duplicates.py --list bgcs_to_remove.txt --folder <gbk_dir>
```

Run the same removal on any FASTA directories, then rebuild the combined FASTA with `gbk2fasta.py` and `concatFasta.py`.

### 7. Map reads (BWA-MEM2)

The reference must end in `.fa` or `.fna`.

```bash
bwa-mem2 index bgcs.fa
sbatch bwa_bam_MT.sh
sbatch bwa-postprocess_bam.sh
```

`bwa-postprocess_bam.sh` runs `bam2txt.py` on each BAM and merges the results into one counts table with `refcode.py`. Reads with equal best hits are split fractionally.

`sam2txt.py` is an alternative counter for SAM input with length, identity and coverage filters (`--min_len`, `--min_identity`, `--min_coverage`).

### 8. Cluster (BiG-SCAPE)

```bash
sbatch bigscape.sh
```

Uses `--alignment-mode auto --include-singletons --gcf-cutoffs 0.6 --mix`.

### 9. Statistics

Fill in the settings at the top of `BGCLandscape_STATS.R` and run it in R. Column order in the counts table must be all case samples, then all controls.

Outputs:
- RPM-normalized counts (wide and long)
- Per-BGC Wilcoxon p-value, rank-biserial correlation (rbc), and permutation FDR (5000 permutations)
- `hit_perm_FDR10`: FDR < 0.10 and |rbc| ≥ 0.3

### 10. Visualize (Cytoscape)

Import the BiG-SCAPE network and use the Organic layout. Filters used:

| Edge attribute | Range |
|---|---|
| distance | 0.3–1.0 |
| DSS score | 0.5–1.0 |

Color nodes by rbc from step 9.

### 11. Taxonomy (MMseqs2)

```bash
python3 extract_contigs_from_bgcs.py \
    --bgc-fasta bgcs.fa \
    --assembly-dir <scaffold_dir> \
    --output-fasta bgc_contigs.fa
sbatch mmseq.sh
```

`mmseq.sh` runs `mmseqs easy-taxonomy` against GTDB with `-s 5 --max-seqs 300 --lca-mode 3 --tax-lineage 1 --orf-filter 0`.
