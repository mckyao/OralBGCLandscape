# RPM normalization, per-BGC Wilcoxon test with rank-biserial correlation,
# and permutation-based empirical FDR.
# Input columns must be ordered: all case samples first, then all controls.

input_file        <- ""   # counts table from bwa-postprocess_bam.sh
mapped_reads_file <- ""   # TSV with columns: sample, mapped_reads
out_norm_wide_csv <- "normalized_RPM_wide.csv"
out_norm_long_csv <- "normalized_RPM_long.csv"
out_stats_csv     <- "bgc_stats_rbc.csv"

n_cavity  <- 27
n_control <- 32

# ---- Load counts and sample metadata ----
mycounts <- read.table(input_file, header = TRUE, row.names = 1, check.names = FALSE)

condition_labels <- c(rep("cavity", n_cavity), rep("control", n_control))
stopifnot(length(condition_labels) == ncol(mycounts))
condition <- factor(condition_labels, levels = c("control", "cavity"))

seq_depth <- read.table(mapped_reads_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("sample", "mapped_reads") %in% colnames(seq_depth)))
seq_depth$sample_base <- sub("_result$", "", seq_depth$sample)

idx <- match(sub("_result$", "", colnames(mycounts)), seq_depth$sample_base)
if (any(is.na(idx))) {
  stop("Samples missing from mapped_reads file:\n  ",
       paste(colnames(mycounts)[is.na(idx)], collapse = ", "))
}

library_sizes <- seq_depth$mapped_reads[idx]
names(library_sizes) <- colnames(mycounts)
if (any(is.na(library_sizes) | library_sizes <= 0)) {
  stop("Samples with missing or zero mapped_reads:\n  ",
       paste(names(library_sizes)[is.na(library_sizes) | library_sizes <= 0], collapse = ", "))
}

# ---- RPM normalization ----
normalized_counts <- sweep(mycounts, 2, library_sizes, "/") * 1e6
write.csv(normalized_counts, out_norm_wide_csv, quote = FALSE)

norm_long <- data.frame(
  feature = rep(rownames(normalized_counts), times = ncol(normalized_counts)),
  sample  = rep(colnames(normalized_counts), each  = nrow(normalized_counts)),
  RPM     = as.vector(as.matrix(normalized_counts)),
  stringsAsFactors = FALSE
)
write.csv(norm_long, out_norm_long_csv, row.names = FALSE, quote = FALSE)

# ---- Wilcoxon test + rank-biserial correlation ----
idx_cavity  <- which(condition == "cavity")
idx_control <- which(condition == "control")
total_pairs <- length(idx_control) * length(idx_cavity)

results_list <- vector("list", nrow(normalized_counts))
for (i in seq_len(nrow(normalized_counts))) {
  counts <- as.numeric(normalized_counts[i, ])
  test <- wilcox.test(counts ~ condition, exact = FALSE)
  U_control <- test$statistic

  # rbc > 0: higher in cavity; rbc < 0: higher in control
  rbc <- 1 - (2 * U_control) / total_pairs

  results_list[[i]] <- data.frame(
    gene           = rownames(normalized_counts)[i],
    U_stat         = U_control,
    p_value        = test$p.value,
    rbc            = rbc,
    mean_caries    = mean(counts[idx_cavity]),
    mean_healthy   = mean(counts[idx_control]),
    median_caries  = median(counts[idx_cavity]),
    median_healthy = median(counts[idx_control]),
    sd_caries      = sd(counts[idx_cavity]),
    sd_healthy     = sd(counts[idx_control])
  )
}
results_df <- do.call(rbind, results_list)

# ---- Permutation-based empirical FDR ----
# mclapply forks on Mac/Linux; on Windows it runs on one core.
set.seed(123)
n_perm     <- 5000
chunk_size <- 200

X <- as.matrix(normalized_counts)
m <- nrow(X)
p_obs <- results_df$p_value

if (anyNA(condition)) stop("condition contains NA")
cond <- as.integer(factor(condition))
if (length(unique(cond)) != 2) stop("condition must have exactly 2 groups")

one_perm <- function(perm_cond, X) {
  g1 <- perm_cond == 1
  g2 <- perm_cond == 2
  vapply(seq_len(nrow(X)), function(i) {
    x1 <- X[i, g1]; x1 <- x1[!is.na(x1)]
    x2 <- X[i, g2]; x2 <- x2[!is.na(x2)]
    if (length(x1) < 1L || length(x2) < 1L) return(NA_real_)
    tryCatch(wilcox.test(x1, x2, exact = FALSE)$p.value, error = function(e) NA_real_)
  }, numeric(1))
}

cores <- max(1L, parallel::detectCores() - 1L)
perm_pvals <- matrix(NA_real_, nrow = m, ncol = n_perm)
pb <- utils::txtProgressBar(min = 0, max = n_perm, style = 3)

perm_ids <- seq_len(n_perm)
for (ids in split(perm_ids, ceiling(perm_ids / chunk_size))) {
  perms <- lapply(ids, function(.) sample(cond))
  res <- parallel::mclapply(perms, function(pc) one_perm(pc, X),
                            mc.cores = cores, mc.set.seed = TRUE)
  mat <- tryCatch(do.call(cbind, res), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) != m || ncol(mat) != length(ids)) {
    stop("Permutation chunk failed")
  }
  perm_pvals[, ids] <- mat
  utils::setTxtProgressBar(pb, max(ids))
}
close(pb)

p_grid <- sort(unique(p_obs))
FDR_empirical <- sapply(p_grid, function(t) {
  fp_null <- mean(colSums(perm_pvals <= t, na.rm = TRUE))
  tp_obs  <- sum(p_obs <= t, na.rm = TRUE)
  if (tp_obs == 0) NA_real_ else fp_null / tp_obs
})

results_df$fdr_perm <- sapply(p_obs, function(p) {
  j <- which(p_grid >= p)[1]
  if (is.na(j)) NA_real_ else FDR_empirical[j]
})

# ---- Hits: empirical FDR < 0.10 and |rbc| >= 0.3 ----
results_df$hit_perm_FDR10 <- !is.na(results_df$fdr_perm) &
  results_df$fdr_perm < 0.10 &
  abs(results_df$rbc) >= 0.3

write.csv(results_df, out_stats_csv, row.names = FALSE)

# ---- Diagnostic plots ----
plot(p_grid, FDR_empirical, type = "l", lwd = 2, col = "blue",
     xlab = "p-value threshold", ylab = "Empirical FDR",
     main = "Permutation-based Empirical FDR")
abline(h = 0.1, col = "red", lty = 2)

perm_flat <- as.vector(perm_pvals)
hist(perm_flat, breaks = 40, col = "grey70",
     main = "Null p-values from permutations", xlab = "p-value")
hist(results_df$p_value, breaks = 40, col = "skyblue",
     main = "Observed p-values", xlab = "p-value")

obs  <- sort(results_df$p_value)
null <- sort(perm_flat[seq_along(obs)])
plot(null, obs, pch = 19, cex = 0.5,
     xlab = "Expected p-values (permutation null)", ylab = "Observed p-values",
     main = "QQ-plot: Observed vs Permutation Null")
abline(0, 1, col = "red")

min_perm <- apply(perm_pvals, 2, min)
hist(min_perm, breaks = 30, col = "grey70",
     main = "Minimum p-value per permutation", xlab = "Minimum permuted p-value")
abline(v = min(results_df$p_value), col = "red", lwd = 2)
legend("topright", legend = "Observed minimum p", col = "red", lwd = 2)
