###############################################################################
# Combinatorial search over haplotype panels for teff panicle shape (1-5 score)
# ---------------------------------------------------------------------------
# The descriptive counterpart of the ordinal regression models, with the same
# two questions as for seed colour: which small sets of haplotypes best
# account for the panicle score, and does the answer depend on the accessions
# used to find it?
#
# PART 1  Best panels on three accession sets. For each panel size k = 1..
#         panel_k_max, every k-haplotype panel is scored on one accession set
#         (the 27 pangenome accessions, the 150 EIAR accessions, or all 177):
#         accessions carrying the same combination of the k haplotypes are
#         given the mean score of that combination, and the panel is scored by
#         the Spearman correlation between this combination mean and the
#         observed score, i.e. how well the panel ranks the accessions from
#         loose to compact. The best panel is recorded, together with the
#         near-best set (all panels within one standard error of the best
#         correlation, SE = 1/sqrt(n - 3)) and the frequency of every
#         haplotype among the near-best panels.
#
# PART 2  Transfer of the 27's best panels to the 150. For each k, the
#         transfer_top_n panels with the highest correlation on the 27 are
#         applied unchanged to the 150: each combination keeps the mean score
#         it had among the 27 (combinations not observed among the 27 receive
#         the overall mean of the 27), and the Spearman correlation between
#         that prediction and the observed score of the 150 is computed.
#
# INPUT (data_file)
#   as for Classifiers_ordinal_panicle.R: Name, code (1-5) and one
#   "haplotype-<locus>" column per locus; haplotypes are written
#   locus:haplotype.
#
# OUTPUTS (under out_dir)
#   COMBINATORIAL_best_panels_by_subset.csv   best panels per k and subset
#   COMBINATORIAL_nearbest_by_subset.csv      haplotype frequencies in the
#                                             near-best sets
#   TRANSFER_27_to_150_panels.csv             every transferred panel
#   TRANSFER_27_to_150_summary.csv            per k: median, quartiles and
#                                             maximum of the correlation on
#                                             the 150
#
# Runtime: panels are scored in vectorised blocks; with about 45 haplotypes
# the search to k = 6 takes of the order of an hour per accession set.
###############################################################################

# =============================================================================
# CONFIGURATION
# =============================================================================
data_file    <- "new_data_panicle.txt"
response_col <- "code"
out_dir      <- "outputs_combinatorial_panicle_ordinal"

# The accession recorded as "DZ" in the data files is DZ_01_354 in the manuscript.
accession_rename <- c(DZ = "DZ_01_354")
train_names_27 <- c(
  "Boni", "DZ_01_354", "Dabbi", "Dtt2-02", "Quncho", "T-33", "T116", "T132",
  "T177", "T206", "T224", "T283", "T288", "T297", "T304", "T330",
  "T336", "T345", "T365", "T366", "T379", "T404", "T412", "T87",
  "T99", "addisie", "karadebi"
)

panel_k_max      <- 6        # largest panel size searched
panel_report_top <- 6        # best panels kept per k and subset
panel_subsets    <- c("pangenome27", "eiar150", "all177")
panel_block_size <- 20000L   # panels scored per block
score_bin        <- 0.001    # resolution at which panel scores are tallied
transfer_top_n   <- 200      # panels of the 27 transferred per k

suppressPackageStartupMessages(library(dplyr))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATA
# =============================================================================
raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
names(raw) <- trimws(names(raw))
stopifnot(all(c("Name", response_col) %in% names(raw)), !anyDuplicated(names(raw)))
raw$Name <- trimws(raw$Name)
hit <- raw$Name %in% names(accession_rename)
raw$Name[hit] <- accession_rename[raw$Name[hit]]
y_all <- suppressWarnings(as.numeric(raw[[response_col]]))
stopifnot(!anyNA(y_all), !anyDuplicated(raw$Name))

hap_raw <- grep("^haplotype[-_]?", names(raw), value = TRUE)
loci    <- sub("^haplotype[-_]?", "", hap_raw)
hap     <- as.data.frame(lapply(raw[hap_raw], as.character), stringsAsFactors = FALSE)
names(hap) <- loci
rownames(hap) <- raw$Name

train_idx <- match(train_names_27, raw$Name)
if (anyNA(train_idx)) stop("Training accessions not found in ", data_file, ": ",
                           paste(train_names_27[is.na(train_idx)], collapse = ", "))
test_idx <- setdiff(seq_len(nrow(hap)), train_idx)
cat("Training accessions:", length(train_idx), "| validation accessions:", length(test_idx),
    "| loci:", length(loci), "\n")

# Binary indicator matrix over the rows `idx`, one column per haplotype at each
# locus, using the haplotypes observed in rows `level_idx`. Columns are
# labelled locus:haplotype.
make_indicators <- function(idx, level_idx = idx) {
  blocks <- lapply(loci, function(m) {
    lv  <- sort(unique(hap[level_idx, m]))
    mat <- matrix(as.integer(outer(hap[idx, m], lv, "==")), ncol = length(lv))
    colnames(mat) <- paste0(m, ":", lv)
    mat
  })
  X <- do.call(cbind, blocks); rownames(X) <- rownames(hap)[idx]
  X
}

# All k-panels starting with column f, as a k x m integer matrix (NULL if none).
panels_from <- function(f, k, P) {
  cands <- if (f < P) (f + 1L):P else integer(0)
  cmb <- if (k == 1L) matrix(seq_len(P), nrow = 1) else
         if (length(cands) < k - 1L) NULL else
         if (length(cands) == k - 1L) matrix(c(f, cands), ncol = 1) else
         rbind(f, combn(cands, k - 1L))
  if (!is.null(cmb)) mode(cmb) <- "integer"
  cmb
}

# =============================================================================
# PART 1 - BEST PANELS ON THE 27, THE 150 AND ALL 177
# =============================================================================
run_descriptive_search <- function(subset_name) {
  idx <- switch(subset_name, pangenome27 = train_idx, eiar150 = test_idx,
                all177 = seq_len(nrow(hap)), stop("unknown subset: ", subset_name))
  X <- make_indicators(idx); y <- y_all[idx]; n <- length(y)
  P <- ncol(X); hap_names <- colnames(X)
  ry <- rank(y); se_rho <- 1 / sqrt(n - 3)
  cat(sprintf("\nPart 1: %s | %d accessions | %d haplotypes | mean score %.2f\n",
              subset_name, n, P, mean(y)))

  codes <- sort(unique(y)); y_code_idx <- match(y, codes)
  nbins <- as.integer(ceiling(1 / score_bin)) + 1L   # correlation in [0, 1] -> bin
  bin_of <- function(s) pmin(nbins, pmax(1L, as.integer(floor(s / score_bin)) + 1L))
  score_of_bin <- (seq_len(nbins) - 1L) * score_bin

  # Spearman correlation of every panel in a block (cmb: k x m column indices)
  # between the combination-mean score and the observed score
  score_block <- function(cmb) {
    k <- nrow(cmb); m <- ncol(cmb); nk <- 2L^k
    K <- matrix(0L, n, m)
    for (j in seq_len(k)) K <- K + X[, cmb[j, ], drop = FALSE] * bitwShiftL(1L, j - 1L)
    offs <- rep(nk * (seq_len(m) - 1L), each = n)
    S <- matrix(0, nk, m); N <- matrix(0L, nk, m)
    for (ci in seq_along(codes)) {
      rows <- y_code_idx == ci
      if (!any(rows)) next
      cnt <- matrix(tabulate(as.vector(K[rows, , drop = FALSE]) + 1L + offs[rep(rows, m)],
                             nbins = nk * m), nrow = nk)
      N <- N + cnt; S <- S + codes[ci] * cnt
    }
    means <- S / pmax(N, 1L)
    pred <- matrix(means[as.vector(K) + 1L + rep(nk * (seq_len(m) - 1L), each = n)], nrow = n)
    rho <- suppressWarnings(as.numeric(cor(apply(pred, 2, rank), ry)))
    rho[is.na(rho)] <- 0
    list(rho = rho, n_pat = colSums(N > 0L))
  }

  top_rows <- list(); near_rows <- list()
  for (k in seq_len(min(panel_k_max, P))) {
    t0 <- Sys.time(); total <- 0
    cnt <- integer(nbins); hmat <- matrix(0L, nbins, P); keep <- NULL
    for (f in if (k == 1L) 1L else seq_len(P - k + 1L)) {
      cmb_f <- panels_from(f, k, P)
      if (is.null(cmb_f)) next
      total <- total + ncol(cmb_f)
      for (b0 in seq.int(1L, ncol(cmb_f), by = panel_block_size)) {
        cmb <- cmb_f[, b0:min(ncol(cmb_f), b0 + panel_block_size - 1L), drop = FALSE]
        sc <- score_block(cmb)
        bn <- bin_of(sc$rho)
        cnt <- cnt + tabulate(bn, nbins = nbins)
        for (j in seq_len(nrow(cmb)))
          hmat <- hmat + matrix(tabulate(bn + nbins * (cmb[j, ] - 1L), nbins = nbins * P), nrow = nbins)
        ord <- head(order(-sc$rho, sc$n_pat), panel_report_top)
        chunk <- data.frame(k = k, spearman = sc$rho[ord], n_patterns = sc$n_pat[ord],
                            accessions_per_pattern = n / sc$n_pat[ord],
                            haplotypes = apply(cmb[, ord, drop = FALSE], 2,
                                               function(z) paste(hap_names[z], collapse = " + ")),
                            stringsAsFactors = FALSE)
        keep <- rbind(keep, chunk)
        keep <- keep[head(order(-keep$spearman, keep$n_patterns), panel_report_top), ]
      }
    }
    best <- max(keep$spearman)
    occupied <- which(cnt > 0L)
    near_cells <- occupied[score_of_bin[occupied] >= best - se_rho - score_bin]
    n_near <- sum(cnt[near_cells]); hcount <- colSums(hmat[near_cells, , drop = FALSE])
    kh <- which(hcount > 0)
    near_rows[[k]] <- data.frame(setting = subset_name, n_accessions = n, k = k, metric = "spearman",
                                 n_total = total, n_near_best = n_near, best = best, tolerance = se_rho,
                                 haplotype = hap_names[kh], n_panels = as.integer(hcount[kh]),
                                 freq = hcount[kh] / n_near, stringsAsFactors = FALSE)
    keep$setting <- subset_name; top_rows[[k]] <- keep
    cat(sprintf("k = %d | %s panels | %.0f s | best rho %.3f | near-best (within 1 SE = %.3f): %s | in >= 50%%: %s\n",
                k, format(total, big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs")),
                best, se_rho, format(n_near, big.mark = ","),
                paste(hap_names[kh][hcount[kh] / n_near >= 0.5], collapse = ", ")))
  }
  list(top = dplyr::bind_rows(top_rows), near = dplyr::bind_rows(near_rows))
}

res <- lapply(panel_subsets, run_descriptive_search)
write.csv(dplyr::bind_rows(lapply(res, `[[`, "top")),
          file.path(out_dir, "COMBINATORIAL_best_panels_by_subset.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(lapply(res, `[[`, "near")),
          file.path(out_dir, "COMBINATORIAL_nearbest_by_subset.csv"), row.names = FALSE)

# =============================================================================
# PART 2 - THE 27's BEST PANELS APPLIED TO THE 150
# =============================================================================
Xtr <- make_indicators(train_idx)
Xte <- make_indicators(test_idx, level_idx = train_idx)   # same columns as Xtr
y_tr <- y_all[train_idx]; y_te <- y_all[test_idx]
ry_tr <- rank(y_tr); ry_te <- rank(y_te)
hap_names <- colnames(Xtr); P <- ncol(Xtr); n_tr <- nrow(Xtr)
cat("\nPart 2: transfer of the 27's top", transfer_top_n, "panels per k to the 150\n")

# Spearman correlation of every panel in a block on the 27
score_train_block <- function(cmb) {
  k <- nrow(cmb); m <- ncol(cmb); nk <- 2L^k
  K <- matrix(0L, n_tr, m)
  for (j in seq_len(k)) K <- K + Xtr[, cmb[j, ], drop = FALSE] * bitwShiftL(1L, j - 1L)
  offs <- nk * (seq_len(m) - 1L)
  N <- matrix(tabulate(as.vector(K) + 1L + rep(offs, each = n_tr), nbins = nk * m), nrow = nk)
  S <- matrix(0, nk, m)
  for (v in sort(unique(y_tr))) {
    rows <- y_tr == v
    S <- S + v * matrix(tabulate(as.vector(K[rows, , drop = FALSE]) + 1L + rep(offs, each = sum(rows)),
                                 nbins = nk * m), nrow = nk)
  }
  means <- S / pmax(N, 1L)
  pred <- matrix(means[as.vector(K) + 1L + rep(offs, each = n_tr)], nrow = n_tr)
  rho <- suppressWarnings(as.numeric(cor(apply(pred, 2, rank), ry_tr))); rho[is.na(rho)] <- 0
  list(rho = rho, n_pat = colSums(N > 0L))
}

rows <- list(); summ <- list()
for (k in seq_len(min(panel_k_max, P))) {
  t0 <- Sys.time(); keep <- NULL
  for (f in if (k == 1L) 1L else seq_len(P - k + 1L)) {
    cmb_f <- panels_from(f, k, P)
    if (is.null(cmb_f)) next
    for (b0 in seq.int(1L, ncol(cmb_f), by = panel_block_size)) {
      cmb <- cmb_f[, b0:min(ncol(cmb_f), b0 + panel_block_size - 1L), drop = FALSE]
      sc <- score_train_block(cmb)
      ord <- head(order(-sc$rho, sc$n_pat), transfer_top_n)
      chunk <- data.frame(k = k, train_spearman = sc$rho[ord], n_patterns_in_27 = sc$n_pat[ord],
                          cols = apply(cmb[, ord, drop = FALSE], 2, function(z) paste(z, collapse = ",")),
                          stringsAsFactors = FALSE)
      keep <- rbind(keep, chunk)
      keep <- keep[head(order(-keep$train_spearman, keep$n_patterns_in_27), transfer_top_n), ]
    }
  }
  # each kept panel: combination means frozen on the 27, correlation on the 150
  keep$test_spearman <- NA_real_; keep$n_test_with_unseen_pattern <- NA_integer_; keep$haplotypes <- NA_character_
  for (i in seq_len(nrow(keep))) {
    cols <- as.integer(strsplit(keep$cols[i], ",")[[1]]); w <- bitwShiftL(1L, seq_along(cols) - 1L)
    ktr <- as.integer(Xtr[, cols, drop = FALSE] %*% w); kte <- as.integer(Xte[, cols, drop = FALSE] %*% w)
    mean_pat <- tapply(y_tr, ktr, mean)
    pred <- unname(mean_pat[as.character(kte)]); unseen <- is.na(pred); pred[unseen] <- mean(y_tr)
    keep$test_spearman[i] <- if (sd(pred) == 0) 0 else suppressWarnings(cor(rank(pred), ry_te))
    keep$n_test_with_unseen_pattern[i] <- sum(unseen)
    keep$haplotypes[i] <- paste(hap_names[cols], collapse = " + ")
  }
  keep$accessions_per_pattern <- n_tr / keep$n_patterns_in_27
  rows[[k]] <- keep[, c("k", "haplotypes", "train_spearman", "n_patterns_in_27", "accessions_per_pattern",
                        "test_spearman", "n_test_with_unseen_pattern")]
  summ[[k]] <- data.frame(k = k, n_panels = nrow(keep),
                          train_spearman_min = min(keep$train_spearman), train_spearman_max = max(keep$train_spearman),
                          test_spearman_median = median(keep$test_spearman),
                          test_spearman_q25 = quantile(keep$test_spearman, 0.25, names = FALSE),
                          test_spearman_q75 = quantile(keep$test_spearman, 0.75, names = FALSE),
                          test_spearman_max = max(keep$test_spearman),
                          share_test_spearman_above_0 = mean(keep$test_spearman > 0))
  cat(sprintf("k = %d | %d panels (rho on the 27: %.2f-%.2f) | rho on the 150: median %+.3f, IQR %+.3f..%+.3f, max %+.3f | %.0f%% > 0 | %.0f s\n",
              k, nrow(keep), min(keep$train_spearman), max(keep$train_spearman), median(keep$test_spearman),
              quantile(keep$test_spearman, .25), quantile(keep$test_spearman, .75), max(keep$test_spearman),
              100 * mean(keep$test_spearman > 0), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
out <- dplyr::bind_rows(rows); out <- out[order(out$k, -out$test_spearman), ]
write.csv(out, file.path(out_dir, "TRANSFER_27_to_150_panels.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(summ), file.path(out_dir, "TRANSFER_27_to_150_summary.csv"), row.names = FALSE)

cat("\nOutputs written under:", normalizePath(out_dir), "\n")
