###############################################################################
# Combinatorial search over haplotype panels for teff seed colour
# ---------------------------------------------------------------------------
# A descriptive counterpart to the classifiers: which small sets of haplotypes
# best account for seed colour, and does the answer depend on the accessions
# used to find it? Two analyses:
#
# PART 1  Perfect separators of the 27 pangenome accessions. Every panel of
#         k = 1..perfect_k_max haplotypes is tested for perfect separation of
#         the 27 into White and Brown (every combination of the k haplotypes
#         carried by accessions of a single colour). Each perfectly separating
#         panel is then applied, unchanged, to the 150 EIAR accessions: the
#         colour assigned to each haplotype combination among the 27 is used to
#         predict the 150, combinations not observed among the 27 receiving
#         the majority colour of the 27, and accuracy is computed on the 150.
#
# PART 2  Best panels on three accession sets. For each panel size k = 1..
#         panel_k_max, every k-haplotype panel is scored on one accession set
#         (the 27, the 150, or all 177) by the number of accessions it assigns
#         correctly when accessions carrying the same combination of the k
#         haplotypes are given the colour prevalent among them. The best panel
#         and the near-best set (all panels within one standard error of the
#         best accuracy) are recorded, together with the frequency of every
#         haplotype among the near-best panels.
#
# INPUT
#   new_data.txt   as for the classifier script: Name, phenotype, one
#                  "haplotype-<locus>" column per locus. Loci 4B/4C are the two
#                  intervals on chromosome 4B (4B1, 4B2), 6A/6B those on
#                  chromosome 6A (6A1, 6A2). Haplotypes are written
#                  locus:haplotype, e.g. 4C:3.
#
# OUTPUTS (under out_dir)
#   PERFECT_SEPARATORS_27_scored_on_150.csv  every perfectly separating panel
#                                            of the 27, with its accuracy on
#                                            the 150
#   PERFECT_SEPARATORS_summary_by_k.csv      count and accuracy range by k
#   COMBINATORIAL_best_panels_by_subset.csv  best panels per k and subset
#   COMBINATORIAL_nearbest_by_subset.csv     haplotype frequencies in the
#                                            near-best sets
#
# Runtime: Part 1 takes minutes. Part 2 is exhaustive; with 69 haplotypes on
# the 177 accessions, k = 5 takes about 15 minutes and k = 6 several hours.
###############################################################################

# =============================================================================
# CONFIGURATION
# =============================================================================
data_file <- "new_data.txt"
out_dir   <- "outputs_combinatorial"

# The accession recorded as "DZ" in new_data.txt is DZ_01_354 in the manuscript.
accession_rename <- c(DZ = "DZ_01_354")

train_names_27 <- c(
  "Boni", "DZ_01_354", "Dabbi", "Dtt2-02", "Quncho", "T-33", "T116", "T132",
  "T177", "T206", "T224", "T283", "T288", "T297", "T304", "T330",
  "T336", "T345", "T365", "T366", "T379", "T404", "T412", "T87",
  "T99", "addisie", "karadebi"
)

perfect_k_max    <- 6      # Part 1: largest panel size searched
panel_k_max      <- 6      # Part 2: largest panel size searched
panel_report_top <- 6      # Part 2: best panels kept per k and subset
panel_subsets    <- c("pangenome27", "eiar150", "all177")

suppressPackageStartupMessages(library(dplyr))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATA
# =============================================================================
raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
names(raw) <- trimws(names(raw))
stopifnot(all(c("Name", "phenotype") %in% names(raw)), !anyDuplicated(names(raw)))
raw$Name <- trimws(raw$Name); raw$phenotype <- trimws(raw$phenotype)
hit <- raw$Name %in% names(accession_rename)
raw$Name[hit] <- accession_rename[raw$Name[hit]]
stopifnot(all(raw$phenotype %in% c("White", "Brown")), !anyDuplicated(raw$Name))

hap_raw <- setdiff(names(raw), c("Name", "phenotype"))
loci    <- sub("^haplotype[-_]?", "", hap_raw)                     # "4B", "4C", ...
hap     <- as.data.frame(lapply(raw[hap_raw], as.character), stringsAsFactors = FALSE)
names(hap) <- loci
rownames(hap) <- raw$Name
pheno <- factor(raw$phenotype, levels = c("White", "Brown"))

train_idx <- match(train_names_27, raw$Name)
if (anyNA(train_idx)) stop("Training accessions not found in ", data_file, ": ",
                           paste(train_names_27[is.na(train_idx)], collapse = ", "))
test_idx <- setdiff(seq_len(nrow(hap)), train_idx)
cat("Training accessions:", length(train_idx), "| validation accessions:", length(test_idx),
    "| loci:", length(loci), "\n")

# Binary indicator matrix: one column per haplotype at each locus, over the
# rows `idx`, using the haplotypes observed in rows `level_idx`. Columns are
# labelled locus:haplotype.
make_indicators <- function(idx, level_idx = idx) {
  blocks <- lapply(loci, function(m) {
    lv  <- sort(unique(hap[level_idx, m]))
    mat <- matrix(as.numeric(outer(hap[idx, m], lv, "==")), ncol = length(lv))
    colnames(mat) <- paste0(m, ":", lv)
    mat
  })
  X <- do.call(cbind, blocks); mode(X) <- "integer"
  rownames(X) <- rownames(hap)[idx]
  X
}

# =============================================================================
# PART 1 - PERFECT SEPARATORS OF THE 27, APPLIED TO THE 150
# =============================================================================
X_train <- make_indicators(train_idx)
X_test  <- make_indicators(test_idx, level_idx = train_idx)   # same columns as X_train
y_train <- pheno[train_idx]; y_test <- pheno[test_idx]
p <- ncol(X_train); haps <- colnames(X_train)
cat("Part 1: haplotypes present among the 27:", p, "\n")

# A panel separates the 27 perfectly if, for every (Brown, White) pair of
# accessions, at least one of its haplotypes differs between the two. Each
# haplotype is encoded as a bitmask over those pairs; a panel is perfect when
# the OR of its masks covers every pair.
B_idx <- which(y_train == "Brown"); W_idx <- which(y_train == "White")
pair_b <- rep(B_idx, each = length(W_idx)); pair_w <- rep(W_idx, times = length(B_idx))
n_pairs <- length(pair_b)
NBITS  <- 31L
nwords <- as.integer(ceiling(n_pairs / NBITS))
word_of <- ((seq_len(n_pairs) - 1L) %/% NBITS) + 1L
bit_of  <- ((seq_len(n_pairs) - 1L) %%  NBITS)
masks <- matrix(0L, nrow = nwords, ncol = p)
for (j in seq_len(p)) {
  diff_pairs <- which(X_train[pair_b, j] != X_train[pair_w, j])
  for (w in unique(word_of[diff_pairs])) {
    bits <- bit_of[diff_pairs[word_of[diff_pairs] == w]]
    masks[w, j] <- Reduce(bitwOr, bitwShiftL(1L, bits), 0L)
  }
}
FULL <- vapply(seq_len(nwords), function(w)
  Reduce(bitwOr, bitwShiftL(1L, bit_of[word_of == w]), 0L), integer(1))
covers <- function(acc) all(acc == FULL)

# All k-panels whose masks cover every pair (running OR over the recursion).
search_perfect_k <- function(k) {
  hits <- list(); nh <- 0L
  rec <- function(start, depth, acc, chosen) {
    if (depth == k - 1L) {
      if (start > p) return(invisible(NULL))
      cand <- start:p
      ok <- rep(TRUE, length(cand))
      for (w in seq_len(nwords)) {
        ok <- ok & (bitwOr(acc[w], masks[w, cand]) == FULL[w])
        if (!any(ok)) return(invisible(NULL))
      }
      for (j in cand[ok]) { nh <<- nh + 1L; hits[[nh]] <<- c(chosen, j) }
      return(invisible(NULL))
    }
    last <- p - (k - depth) + 1L
    if (start > last) return(invisible(NULL))
    for (j in start:last) rec(j + 1L, depth + 1L, bitwOr(acc, masks[, j]), c(chosen, j))
    invisible(NULL)
  }
  rec(1L, 0L, integer(nwords), integer(0))
  hits
}

# Minimal: no haplotype can be dropped without losing the separation.
is_minimal <- function(cols) {
  for (d in seq_along(cols)) {
    acc <- integer(nwords)
    for (j in cols[-d]) acc <- bitwOr(acc, masks[, j])
    if (covers(acc)) return(FALSE)
  }
  TRUE
}

# Apply a panel found on the 27 to the 150: each combination of the panel's
# haplotypes takes the majority colour of the 27 accessions carrying it
# (ties and combinations absent from the 27 take the majority colour of the 27).
train_majority   <- names(which.max(table(y_train)))
fallback_is_brown <- train_majority == "Brown"
brown_tr <- y_train == "Brown"; brown_te <- y_test == "Brown"
nB_test <- sum(brown_te); nW_test <- sum(!brown_te)

score_on_test <- function(cols) {
  k <- length(cols); w <- as.integer(2^(seq_len(k) - 1L)); nk <- 2L^k
  ktr <- as.integer(X_train[, cols, drop = FALSE] %*% w) + 1L
  kte <- as.integer(X_test[,  cols, drop = FALSE] %*% w) + 1L
  nb <- tabulate(ktr[brown_tr], nbins = nk); nw <- tabulate(ktr[!brown_tr], nbins = nk)
  seen_pat  <- (nb + nw) > 0L
  pat_brown <- ifelse(nb > nw, TRUE, ifelse(nw > nb, FALSE, fallback_is_brown))
  pat_brown[!seen_pat] <- fallback_is_brown
  pred_brown <- pat_brown[kte]
  TP <- sum(brown_te & pred_brown); TN <- sum(!brown_te & !pred_brown)
  c(acc = (TP + TN) / length(kte), bal = mean(c(TP / nB_test, TN / nW_test)),
    sens = TP / nB_test, spec = TN / nW_test,
    n_unseen = sum(!seen_pat[kte]), n_patterns = sum(seen_pat))
}

perfect <- list()
for (k in seq_len(min(perfect_k_max, p))) {
  t0 <- Sys.time()
  hits <- search_perfect_k(k)
  cat(sprintf("k = %d | %s panels | %d perfectly separate the 27 (%.0f s)\n",
              k, format(choose(p, k), big.mark = ",", scientific = FALSE), length(hits),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (!length(hits)) next
  perfect[[k]] <- dplyr::bind_rows(lapply(hits, function(h) {
    sc <- score_on_test(h)
    data.frame(k = k, minimal = is_minimal(h),
               haplotypes = paste(haps[h], collapse = " + "),
               n_patterns_in_27 = sc[["n_patterns"]],
               accessions_per_pattern = length(train_idx) / sc[["n_patterns"]],
               test_accuracy = sc[["acc"]], test_balanced_accuracy = sc[["bal"]],
               test_sensitivity_brown = sc[["sens"]], test_specificity_white = sc[["spec"]],
               n_test_with_unseen_pattern = sc[["n_unseen"]],
               stringsAsFactors = FALSE)
  }))
}
if (!length(perfect)) stop("No perfectly separating panel of up to ", perfect_k_max, " haplotypes.")
perfect <- dplyr::bind_rows(perfect)
perfect <- perfect[order(perfect$k, -perfect$test_accuracy), ]
write.csv(perfect, file.path(out_dir, "PERFECT_SEPARATORS_27_scored_on_150.csv"), row.names = FALSE)

summary_k <- perfect %>%
  dplyr::group_by(k, minimal) %>%
  dplyr::summarise(n_panels = dplyr::n(),
                   min_test_accuracy = min(test_accuracy),
                   median_test_accuracy = median(test_accuracy),
                   max_test_accuracy = max(test_accuracy), .groups = "drop")
write.csv(summary_k, file.path(out_dir, "PERFECT_SEPARATORS_summary_by_k.csv"), row.names = FALSE)

k_min <- min(perfect$k)
cat(sprintf("\nMinimum panel size separating the 27 perfectly: k = %d (%d panels)\n",
            k_min, sum(perfect$k == k_min)))
print(perfect[perfect$k == k_min, c("haplotypes", "n_patterns_in_27", "test_accuracy",
                                    "n_test_with_unseen_pattern")],
      row.names = FALSE, digits = 3)
cat("\nBy panel size:\n"); print(as.data.frame(summary_k), row.names = FALSE, digits = 3)

# =============================================================================
# PART 2 - BEST PANELS ON THE 27, THE 150 AND ALL 177
# =============================================================================
run_descriptive_search <- function(subset_name) {
  idx <- switch(subset_name,
                pangenome27 = train_idx, eiar150 = test_idx,
                all177 = seq_len(nrow(hap)), stop("unknown subset: ", subset_name))
  X <- make_indicators(idx)               # haplotypes present in this subset
  n <- nrow(X); P <- ncol(X); hap_names <- colnames(X)
  y <- as.integer(pheno[idx] == "Brown"); yw <- 1L - y
  cat(sprintf("\nPart 2: %s | %d accessions | %d haplotypes | majority-rule baseline %.3f\n",
              subset_name, n, P, max(mean(y), 1 - mean(y))))

  # accessions correctly assigned when every combination of the panel's
  # haplotypes takes the colour prevalent among its carriers
  n_correct <- function(cols) {
    k <- length(cols); w <- as.integer(2^(seq_len(k) - 1L)); nk <- 2L^k
    key <- as.integer(X[, cols, drop = FALSE] %*% w) + 1L
    nb <- tabulate(key[y == 1L], nbins = nk); nw <- tabulate(key[yw == 1L], nbins = nk)
    c(correct = sum(pmax(nb, nw)), n_pat = sum((nb + nw) > 0L))
  }

  top_rows <- list(); near_rows <- list()
  for (k in seq_len(panel_k_max)) {
    t0 <- Sys.time(); total <- 0
    cnt <- integer(n + 1L); hmat <- matrix(0L, nrow = n + 1L, ncol = P)
    best <- -1L; keep <- NULL
    firsts <- if (k == 1L) 1L else seq_len(P - k + 1L)
    for (f in firsts) {
      cands <- if (f < P) (f + 1L):P else integer(0)
      cmb <- if (k == 1L) matrix(seq_len(P), nrow = 1) else
             if (length(cands) == 1L) matrix(c(f, cands), ncol = 1) else
             rbind(f, combn(cands, k - 1L))
      res <- apply(cmb, 2, n_correct); total <- total + ncol(cmb)
      corr <- res["correct", ]
      best <- max(best, max(corr))
      # tally panels by score (4 accessions is the widest tolerance used below)
      hit <- which(corr >= best - 4L)
      for (u in unique(corr[hit])) {
        rows_u <- hit[corr[hit] == u]
        cnt[u + 1L] <- cnt[u + 1L] + length(rows_u)
        hmat[u + 1L, ] <- hmat[u + 1L, ] + tabulate(as.vector(cmb[, rows_u, drop = FALSE]), nbins = P)
      }
      ord <- head(order(-corr, res["n_pat", ]), panel_report_top)
      chunk <- data.frame(k = k, accuracy = corr[ord] / n, n_patterns = res["n_pat", ord],
                          accessions_per_pattern = n / res["n_pat", ord],
                          haplotypes = apply(cmb[, ord, drop = FALSE], 2,
                                             function(z) paste(hap_names[z], collapse = " + ")),
                          stringsAsFactors = FALSE)
      keep <- rbind(keep, chunk)
      keep <- keep[head(order(-keep$accuracy, keep$n_patterns), panel_report_top), ]
    }
    # near-best set: within one standard error of the best accuracy, expressed
    # in accessions (at least one)
    pb  <- best / n
    tol <- max(1L, as.integer(round(sqrt(pb * (1 - pb) / n) * n)))
    lv  <- (max(0L, best - tol):best) + 1L
    n_near <- sum(cnt[lv]); hcount <- colSums(hmat[lv, , drop = FALSE])
    kh <- which(hcount > 0)
    near_rows[[k]] <- data.frame(setting = subset_name, n_accessions = n, k = k,
                                 n_total = total, n_near_best = n_near, best = pb,
                                 tolerance_accessions = tol,
                                 haplotype = hap_names[kh], n_panels = as.integer(hcount[kh]),
                                 freq = hcount[kh] / n_near, stringsAsFactors = FALSE)
    keep$setting <- subset_name; top_rows[[k]] <- keep
    cat(sprintf("k = %d | %s panels | %.0f s | best %.3f | near-best (within %d accessions): %s | in >= 50%%: %s\n",
                k, format(total, big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "secs")), pb, tol,
                format(n_near, big.mark = ","),
                paste(hap_names[kh][hcount[kh] / n_near >= 0.5], collapse = ", ")))
  }
  list(top = dplyr::bind_rows(top_rows), near = dplyr::bind_rows(near_rows))
}

res <- lapply(panel_subsets, run_descriptive_search)
write.csv(dplyr::bind_rows(lapply(res, `[[`, "top")),
          file.path(out_dir, "COMBINATORIAL_best_panels_by_subset.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(lapply(res, `[[`, "near")),
          file.path(out_dir, "COMBINATORIAL_nearbest_by_subset.csv"), row.names = FALSE)

cat("\nOutputs written under:", normalizePath(out_dir), "\n")
