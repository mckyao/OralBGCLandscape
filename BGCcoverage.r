# Flags BGCs that are >=95% covered by a longer BGC in the all-vs-all BLAST.
# Input: *_ambiguous_windows.*.tsv from run_bgc_blast.sh

library(dplyr)
library(stringr)

input_file <- ""   # ambiguous windows TSV
cov_cutoff <- 0.95

df <- read.table(input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("qseqid", "sseqid", "pident", "aln_len", "qlen",
                               "qstart", "qend", "sstart", "send")) %>%
  mutate(
    qstart    = as.numeric(qstart),
    qend      = as.numeric(qend),
    qlen_true = as.numeric(str_extract(qseqid, "(?<=length_)\\d+")),
    slen_true = as.numeric(str_extract(sseqid, "(?<=length_)\\d+"))
  )

merge_intervals <- function(starts, ends) {
  ok <- !(is.na(starts) | is.na(ends))
  starts <- as.numeric(starts[ok]); ends <- as.numeric(ends[ok])
  if (length(starts) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
  ord <- order(starts, ends)
  starts <- starts[ord]; ends <- ends[ord]
  ms <- c(); me <- c()
  cs <- starts[1]; ce <- ends[1]
  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= ce) {
        ce <- max(ce, ends[i])
      } else {
        ms <- c(ms, cs); me <- c(me, ce)
        cs <- starts[i]; ce <- ends[i]
      }
    }
  }
  data.frame(start = c(ms, cs), end = c(me, ce))
}

covered_bp <- function(s, e) {
  iv <- merge_intervals(s, e)
  if (nrow(iv) == 0) 0 else sum(iv$end - iv$start + 1)
}

# ---- Per-BGC coverage summary ----
bgc_overlap_table <- df %>%
  filter(!is.na(qstart), !is.na(qend), !is.na(qlen_true)) %>%
  group_by(qseqid, qlen_true) %>%
  summarise(ambiguous_bp = covered_bp(qstart, qend), .groups = "drop") %>%
  mutate(percent_covered = 100 * ambiguous_bp / qlen_true) %>%
  rename(bgc = qseqid) %>%
  arrange(desc(percent_covered))

amb <- bgc_overlap_table$percent_covered[bgc_overlap_table$percent_covered > 0]
cat("Total BGCs:", nrow(bgc_overlap_table), "\n")
cat("BGCs with any overlap:", length(amb), "\n")
if (length(amb) > 0) {
  cat("Median % covered:", signif(median(amb), 3), "\n")
  cat("IQR % covered:", signif(quantile(amb, c(0.25, 0.75)), 3), "\n")
  cat("Max % covered:", signif(max(amb), 3), "\n")
}
cat("BGCs >=95% covered:", sum(bgc_overlap_table$percent_covered >= 95), "\n")

write.csv(bgc_overlap_table, "bgc_overlap_percent_table.csv", row.names = FALSE)

# ---- Shorter near-duplicates to remove ----
to_remove <- df %>%
  filter(!is.na(qstart), !is.na(qend), !is.na(qlen_true), !is.na(slen_true)) %>%
  group_by(qseqid, sseqid, qlen_true, slen_true) %>%
  summarise(ambiguous_bp = covered_bp(qstart, qend), .groups = "drop") %>%
  filter(ambiguous_bp / qlen_true >= cov_cutoff, qlen_true < slen_true) %>%
  pull(qseqid) %>%
  unique()

cat("BGCs to remove:", length(to_remove), "\n")
writeLines(to_remove, "bgcs_to_remove.txt")
