###############################################################################
# Haplotype-based prediction of teff panicle shape as an ordinal trait
# ---------------------------------------------------------------------------
# Panicle form is scored in the field on a five-point scale (1 very loose ...
# 5 very compact). The score is modelled as a number by twelve regression
# variants trained on the 27 pangenome accessions and evaluated (i) internally
# on the 27, (ii) externally on the 150 EIAR accessions never used for
# training or tuning, and (iii) by 5-fold cross-validation over all 177
# accessions. Haplotype importance is summarised as a rank consensus across
# the variants whose external Spearman correlation exceeds a threshold, and
# every haplotype is tested for association with the score on all 177
# accessions.
#
# INPUT (data_file)
#   tab-separated, one row per accession: Name, code (1-5) and one
#   "haplotype-<locus>" column per locus, with haplotypes called jointly on all
#   177 accessions. The locus id is the column name without the prefix
#   (<chromosome>_<window start>_<window end>). Other columns are ignored.
#
# MODELS (each fitted with and without predictor standardisation)
#   L1- and L2-regularised linear regression (glmnet, gaussian), linear and
#   polynomial-kernel support vector regression (e1071, eps-regression),
#   regression trees (rpart, anova) and gradient boosting (gbm, gaussian).
#   Haplotypes enter as binary indicators (one per haplotype at each locus).
#   Hyperparameters are tuned by inner cross-validation on the training rows
#   of each fold, with folds balanced over the score and the Spearman
#   correlation between predicted and observed score as the criterion;
#   standardisation uses training-fold statistics only; indicators constant in
#   the training rows of a fold are removed before fitting.
#
# EVALUATION
#   Spearman correlation between predicted and observed score (primary; the
#   correlation of a constant prediction is scored as 0), with its P value,
#   plus MAE and RMSE in score units and the MAE of the constant predictor
#   (median score of the training fold).
#   Regimes: resub27 (fitted and scored on the 27), loocv27 (leave-one-out on
#   the 27), cv5_27 (5-fold on the 27), holdout150 (fitted on the 27, scored
#   on the 150), cv5_all (5-fold over all 177). Each regime is repeated for
#   every seed in seeds_to_run and averaged.
#
# UNSEEN HAPLOTYPES
#   A haplotype present among the 150 but absent from the 27 has no indicator
#   in a model fitted on the 27; accessions carrying it are kept and flagged,
#   and the holdout150 correlation is also reported on the accessions
#   carrying only training-set haplotypes.
#
# ASSOCIATION (no model, all 177 accessions)
#   For every haplotype, the score of its carriers is compared with that of
#   all other accessions by the Mann-Whitney test, and the share of very
#   compact panicles (code >= very_compact_from) among carriers with all
#   others by Fisher's exact test; both are Benjamini-Hochberg adjusted.
#
# OUTPUTS (under out_dir)
#   accessions_split.csv, unseen_haplotypes.csv      as in the binary script
#   <regime>/seed_<s>/fold_results.csv               per-accession predictions
#   <regime>/seed_<s>/<model>_importance.csv         per-haplotype importance
#   <regime>/mean_metrics_over_seeds.csv             metrics averaged over seeds
#   METRICS_by_model_and_regime.csv                  all regimes
#   METRICS_spearman_wide.csv                        model x regime
#   IMPORTANCE_ranks_long_holdout150.csv             rank of every haplotype in
#                                                    every retained model
#   IMPORTANCE_consensus_holdout150.csv              median rank, number of
#                                                    models, direction
#   HAPLOTYPE_association.csv                        Mann-Whitney and
#                                                    very-compact contrast
#
# Runtime: minutes per seed for the 27-accession regimes; cv5_all is the
# slowest (training sets of ~140 rows).
###############################################################################

# =============================================================================
# CONFIGURATION
# =============================================================================
data_file    <- "new_data_panicle.txt"
response_col <- "code"
out_dir      <- "outputs_classifiers_panicle_ordinal"
dir_high     <- "Compact"     # label for a positive effect (higher score)
dir_low      <- "Loose"       # label for a negative effect (lower score)
very_compact_from <- 4        # association contrast: score >= this vs the rest

# The accession recorded as "DZ" in the data files is DZ_01_354 in the manuscript.
accession_rename <- c(DZ = "DZ_01_354")
train_names_27 <- c(
  "Boni", "DZ_01_354", "Dabbi", "Dtt2-02", "Quncho", "T-33", "T116", "T132",
  "T177", "T206", "T224", "T283", "T288", "T297", "T304", "T330",
  "T336", "T345", "T365", "T366", "T379", "T404", "T412", "T87",
  "T99", "addisie", "karadebi"
)

seeds_to_run   <- 1:3
regimes_to_run <- c("resub27", "loocv27", "cv5_27", "holdout150", "cv5_all")
cv_folds       <- 5
innerK_base    <- 5

# Hyperparameter grids
svr_cost_grid        <- c(0.01, 0.1, 1, 10)
svr_epsilon          <- 0.1
svr_poly_cost_grid   <- c(0.01, 0.1, 1)
svr_poly_degree_grid <- c(2, 3)
svr_poly_coef0       <- 1
cart_cp_grid         <- c(0.0, 0.001, 0.01, 0.05)
cart_maxdepth        <- 30
cart_minsplit        <- 2
gbm_depth_grid       <- c(1, 2, 3)
gbm_trees_grid       <- c(50, 100, 200, 500)
gbm_shrinkage        <- 0.01
gbm_minobs           <- 2
gbm_bagfrac          <- 0.8
svr_poly_nperm       <- 10    # permutations for the polynomial-SVR importance

# Importance consensus: models retained when their mean external Spearman
# correlation exceeds the threshold (0.16 is two standard errors above zero
# for 150 accessions, SE = 1/sqrt(n - 3) = 0.082).
importance_regime    <- "holdout150"
importance_threshold <- 0.16

# Association: haplotypes with fewer carriers (or fewer non-carriers) are
# listed but not tested.
assoc_min_carriers <- 3
assoc_alpha        <- 0.05

model_names <- c("lasso_nostd", "lasso_std", "ridge_nostd", "ridge_std",
                 "svr_lin_nostd", "svr_lin_std", "svr_poly_nostd", "svr_poly_std",
                 "cart_nostd", "cart_std", "gbm_nostd", "gbm_std")
signed_models <- c("lasso_nostd", "lasso_std", "ridge_nostd", "ridge_std",
                   "svr_lin_nostd", "svr_lin_std")

pkgs <- c("dplyr", "tidyr", "glmnet", "e1071", "rpart", "gbm")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# METRICS AND HELPERS
# =============================================================================
safe_cor <- function(a, b) {
  if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(0)
  r <- suppressWarnings(cor(a, b, method = "spearman"))
  if (is.na(r)) 0 else r
}
safe_cor_p <- function(a, b) {
  if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)$p.value)
}

# Standardise with training statistics only.
prep_std <- function(Xtr, Xva, standardize_flag) {
  if (!standardize_flag) return(list(Xtr = Xtr, Xva = Xva, sigma = rep(1, ncol(Xtr))))
  mu <- colMeans(Xtr); sigma <- apply(Xtr, 2, sd)
  sigma[is.na(sigma) | sigma == 0] <- 1
  list(Xtr = sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/"),
       Xva = sweep(sweep(Xva, 2, mu, "-"), 2, sigma, "/"), sigma = sigma)
}

# Folds balanced over the score: accessions are sorted by score and dealt out
# in turn, so every fold sees the rare scores whenever there are enough of them.
make_folds_ordinal <- function(y, K = 5, seed = 1) {
  set.seed(seed)
  n <- length(y)
  o <- order(y, runif(n))
  fold <- integer(n)
  fold[o] <- rep(sample(K), length.out = n)
  fold
}

# =============================================================================
# INNER-CV TUNERS (criterion: Spearman correlation)
# =============================================================================
tune_glmnet_innerK <- function(Xtr, ytr, K, standardize_flag, alpha_glmnet, seed = 1) {
  folds <- make_folds_ordinal(ytr, K, seed)
  fit0 <- glmnet::glmnet(x = Xtr, y = ytr, family = "gaussian", alpha = alpha_glmnet,
                         standardize = standardize_flag)
  lambda_grid <- fit0$lambda
  fold_scores <- matrix(NA_real_, K, length(lambda_grid))
  for (k in 1:K) {
    va <- which(folds == k); tr <- which(folds != k)
    if (length(va) < 2) next
    fit <- glmnet::glmnet(x = Xtr[tr, , drop = FALSE], y = ytr[tr], family = "gaussian",
                          alpha = alpha_glmnet, lambda = lambda_grid, standardize = standardize_flag)
    pmat <- predict(fit, newx = Xtr[va, , drop = FALSE])
    idx_used <- match(fit$lambda, lambda_grid)
    for (j in seq_along(fit$lambda))
      fold_scores[k, idx_used[j]] <- safe_cor(as.numeric(pmat[, j]), ytr[va])
  }
  # one-standard-error rule: the most regularised lambda within one SE of the best
  mean_sc <- colMeans(fold_scores, na.rm = TRUE)
  sd_sc   <- apply(fold_scores, 2, sd, na.rm = TRUE)
  se_sc   <- sd_sc / sqrt(pmax(colSums(!is.na(fold_scores)), 1)); se_sc[is.na(se_sc)] <- 0
  j_best  <- which.max(mean_sc)
  cand    <- which(mean_sc >= mean_sc[j_best] - se_sc[j_best])
  if (!length(cand)) cand <- j_best
  list(lambda = lambda_grid[max(cand)])
}

fit_svr <- function(X, y, kernel, cost, degree = NULL) {
  if (kernel == "linear") {
    e1071::svm(x = X, y = y, type = "eps-regression", kernel = "linear",
               cost = cost, epsilon = svr_epsilon, scale = FALSE)
  } else {
    e1071::svm(x = X, y = y, type = "eps-regression", kernel = "polynomial",
               cost = cost, degree = degree, coef0 = svr_poly_coef0,
               epsilon = svr_epsilon, scale = FALSE)
  }
}

tune_svr_innerK <- function(Xtr, ytr, K, standardize_flag, cost_grid, kernel,
                            degree_grid = NULL, seed = 1) {
  folds <- make_folds_ordinal(ytr, K, seed)
  grid <- if (kernel == "linear") data.frame(cost = cost_grid, degree = NA) else
    expand.grid(cost = cost_grid, degree = degree_grid)
  scores <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(va) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_svr(st$Xtr, ytr[tr], kernel, grid$cost[g], grid$degree[g])
      fold_sc[k] <- safe_cor(as.numeric(predict(fit, st$Xva)), ytr[va])
    }
    scores[g] <- mean(fold_sc, na.rm = TRUE)
  }
  b <- which.max(scores)
  list(cost = grid$cost[b], degree = grid$degree[b])
}

fit_cart <- function(X, y, cp) {
  df_tr <- data.frame(y = y, as.data.frame(X))
  rpart::rpart(y ~ ., data = df_tr, method = "anova",
               control = rpart::rpart.control(cp = cp, maxdepth = cart_maxdepth,
                                              minsplit = cart_minsplit, xval = 0))
}

tune_cart_innerK <- function(Xtr, ytr, K, standardize_flag, cp_grid, seed = 1) {
  folds <- make_folds_ordinal(ytr, K, seed)
  scores <- setNames(rep(NA_real_, length(cp_grid)), as.character(cp_grid))
  for (cpv in cp_grid) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(va) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_cart(st$Xtr, ytr[tr], cpv)
      fold_sc[k] <- safe_cor(as.numeric(predict(fit, newdata = as.data.frame(st$Xva))), ytr[va])
    }
    scores[as.character(cpv)] <- mean(fold_sc, na.rm = TRUE)
  }
  list(cp = as.numeric(names(which.max(scores))))
}

fit_gbm <- function(X, y, depth, seed) {
  df_tr <- data.frame(y = y, as.data.frame(X))
  set.seed(seed)
  gbm::gbm(y ~ ., data = df_tr, distribution = "gaussian", n.trees = max(gbm_trees_grid),
           interaction.depth = depth, shrinkage = gbm_shrinkage, n.minobsinnode = gbm_minobs,
           bag.fraction = gbm_bagfrac, train.fraction = 1.0, cv.folds = 0, verbose = FALSE)
}

tune_gbm_innerK <- function(Xtr, ytr, K, standardize_flag, depth_grid, trees_grid, seed = 1) {
  folds <- make_folds_ordinal(ytr, K, seed)
  grid <- expand.grid(depth = depth_grid, trees = trees_grid)
  scores <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(va) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_gbm(st$Xtr, ytr[tr], grid$depth[g],
                     seed = seed + 10000 * g + 100 * k + if (standardize_flag) 1 else 0)
      p_va <- predict(fit, newdata = as.data.frame(st$Xva), n.trees = grid$trees[g])
      fold_sc[k] <- safe_cor(as.numeric(p_va), ytr[va])
    }
    scores[g] <- mean(fold_sc, na.rm = TRUE)
  }
  b <- which.max(scores)
  list(depth = grid$depth[b], trees = grid$trees[b])
}

# =============================================================================
# IMPORTANCE EXTRACTORS
# =============================================================================
# Signed (positive = higher score): coefficients of the regularised linear
# models and weights of the linear SVR. Unsigned: permutation drop in Spearman
# correlation (polynomial SVR), variable importance (CART), relative influence
# (GBM).
importance_glmnet <- function(fit, lambda, feat_names) {
  b <- as.matrix(coef(fit, s = lambda))
  b <- b[setdiff(rownames(b), "(Intercept)"), , drop = FALSE]
  v <- setNames(rep(0, length(feat_names)), feat_names)
  common <- intersect(rownames(b), feat_names)
  v[common] <- b[common, 1]
  v
}

importance_svr_linear <- function(svm_fit, Xtr_fit, sigma) {
  w <- drop(t(svm_fit$coefs) %*% svm_fit$SV)
  names(w) <- colnames(Xtr_fit)
  w / sigma                       # back to the unstandardised scale
}

importance_perm_rho <- function(fit, Xtr_fit, ytr, nperm = 10, seed = 1) {
  set.seed(seed)
  base <- safe_cor(as.numeric(predict(fit, Xtr_fit)), ytr)
  imp <- setNames(rep(0, ncol(Xtr_fit)), colnames(Xtr_fit))
  for (j in seq_len(ncol(Xtr_fit))) {
    drops <- numeric(nperm)
    for (b in 1:nperm) {
      Xp <- Xtr_fit; Xp[, j] <- sample(Xp[, j])
      drops[b] <- base - safe_cor(as.numeric(predict(fit, Xp)), ytr)
    }
    imp[j] <- mean(pmax(drops, 0))
  }
  imp
}

importance_cart_vi <- function(fit, feat_names) {
  v <- setNames(rep(0, length(feat_names)), feat_names)
  vi <- fit$variable.importance
  if (!is.null(vi) && length(vi)) {
    common <- intersect(names(vi), feat_names); v[common] <- as.numeric(vi[common])
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

importance_gbm_rel <- function(fit, n.trees, feat_names) {
  v <- setNames(rep(0, length(feat_names)), feat_names)
  rel <- suppressWarnings(summary(fit, plotit = FALSE, n.trees = n.trees))
  if (!is.null(rel) && nrow(rel)) {
    common <- intersect(as.character(rel$var), feat_names)
    v[common] <- rel$rel.inf[match(common, as.character(rel$var))]
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

# Fold-wise importances of one model -> one row per haplotype. `magnitude` is
# the absolute mean across folds; `direction` (signed models) is its sign.
summarise_importance <- function(M, signed, label_of) {
  mean_i <- colMeans(M, na.rm = TRUE)
  out <- data.frame(feature = colnames(M), haplotype = unname(label_of[colnames(M)]),
                    mean = as.numeric(mean_i), sd = as.numeric(apply(M, 2, sd, na.rm = TRUE)),
                    magnitude = abs(as.numeric(mean_i)), stringsAsFactors = FALSE)
  if (signed) out$direction <- ifelse(mean_i > 0, dir_high, dir_low)
  out[order(-out$magnitude), ]
}

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
stopifnot(!anyNA(y_all), all(y_all >= 1 & y_all <= 5), !anyDuplicated(raw$Name))

hap_raw    <- grep("^haplotype[-_]?", names(raw), value = TRUE)
haplo_cols <- sub("^haplotype[-_]?", "", hap_raw)                  # locus id
df0 <- data.frame(row.names = raw$Name)
for (j in seq_along(hap_raw)) df0[[haplo_cols[j]]] <- factor(as.character(raw[[hap_raw[j]]]))

train_idx <- match(train_names_27, raw$Name)
if (anyNA(train_idx)) stop("Training accessions not found in ", data_file, ": ",
                           paste(train_names_27[is.na(train_idx)], collapse = ", "))
test_idx <- setdiff(seq_len(nrow(df0)), train_idx)

# One binary indicator per haplotype at each locus, over all accessions.
# Column names are syntactically valid R names (needed by rpart and gbm);
# `feat_label` carries the display form locus:haplotype.
X_all <- do.call(cbind, lapply(haplo_cols, function(m) {
  lv  <- levels(df0[[m]])
  mat <- matrix(as.numeric(outer(as.character(df0[[m]]), lv, "==")), ncol = length(lv))
  colnames(mat) <- make.names(paste0(m, "_", lv))
  mat
}))
rownames(X_all) <- rownames(df0)
feat_names   <- colnames(X_all)
feat_label   <- unlist(lapply(haplo_cols, function(m) paste0(m, ":", levels(df0[[m]]))))
label_of     <- setNames(feat_label, feat_names)
locus_of_col <- rep(haplo_cols, vapply(haplo_cols, function(m) nlevels(df0[[m]]), integer(1)))
ids <- rownames(df0)

present_in_27 <- colSums(X_all[train_idx, , drop = FALSE]) > 0
unseen_mat <- vapply(haplo_cols, function(m) {
  cols <- which(locus_of_col == m & present_in_27)
  rowSums(X_all[test_idx, cols, drop = FALSE]) == 0
}, logical(length(test_idx)))
has_unseen_level <- rep(NA, nrow(df0))
has_unseen_level[test_idx] <- rowSums(unseen_mat) > 0

cat("Training accessions:", length(train_idx), "| validation accessions:", length(test_idx),
    "| loci:", length(haplo_cols), "| haplotype indicators:", ncol(X_all), "\n")
cat("Score distribution, 27: "); print(table(y_all[train_idx]))
cat("Score distribution, 150:"); print(table(y_all[test_idx]))

write.csv(data.frame(accession = ids,
                     set = ifelse(seq_along(ids) %in% train_idx, "train27", "test150"),
                     code = y_all, has_unseen_haplotype = has_unseen_level),
          file.path(out_dir, "accessions_split.csv"), row.names = FALSE)
write.csv(data.frame(locus = haplo_cols,
                     haplotypes_in_27 = vapply(haplo_cols, function(m)
                       paste(levels(droplevels(df0[[m]][train_idx])), collapse = "/"), ""),
                     n_test_with_unseen_haplotype = as.integer(colSums(unseen_mat))),
          file.path(out_dir, "unseen_haplotypes.csv"), row.names = FALSE)

# =============================================================================
# FOLD PLANS
# =============================================================================
make_fold_plan <- function(regime, seed) {
  tr <- train_idx
  switch(regime,
    resub27    = list(list(train = tr, test = tr)),
    loocv27    = lapply(seq_along(tr), function(j) list(train = tr[-j], test = tr[j])),
    cv5_27     = {
      fid <- make_folds_ordinal(y_all[tr], cv_folds, seed)
      lapply(sort(unique(fid)), function(k) list(train = tr[fid != k], test = tr[fid == k]))
    },
    holdout150 = list(list(train = tr, test = test_idx)),
    cv5_all    = {
      all_idx <- seq_along(y_all); fid <- make_folds_ordinal(y_all, cv_folds, seed)
      lapply(sort(unique(fid)), function(k) list(train = all_idx[fid != k], test = all_idx[fid == k]))
    },
    stop("unknown regime: ", regime))
}

# =============================================================================
# ONE RUN: ONE SEED, ONE REGIME
# =============================================================================
run_one_seed <- function(seed, regime, run_dir) {
  set.seed(seed)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  n <- nrow(X_all)
  fold_plan  <- make_fold_plan(regime, seed)
  n_folds    <- length(fold_plan)
  tested_idx <- sort(unique(unlist(lapply(fold_plan, function(f) f$test))))

  pred_store <- setNames(lapply(model_names, function(x) rep(NA_real_, n)), model_names)
  IMP <- setNames(lapply(model_names, function(x)
    matrix(0, n_folds, length(feat_names),
           dimnames = list(paste0("fold", seq_len(n_folds)), feat_names))), model_names)
  train_median <- rep(NA_real_, n)     # per tested row: median score of its training fold
  store_imp <- function(tag, i, v_keep, feat_outer) {
    v_full <- setNames(rep(0, length(feat_names)), feat_names)
    v_full[feat_outer] <- v_keep
    IMP[[tag]][i, ] <<- v_full
  }

  for (i in seq_len(n_folds)) {
    te <- fold_plan[[i]]$test; tr <- fold_plan[[i]]$train
    ytr <- y_all[tr]
    train_median[te] <- median(ytr)
    innerK <- max(2, min(innerK_base, length(ytr) %/% 3))

    sd_outer <- apply(X_all[tr, , drop = FALSE], 2, sd); sd_outer[is.na(sd_outer)] <- 0
    keep_outer <- sd_outer > 0
    if (!any(keep_outer)) {
      for (nm in model_names) pred_store[[nm]][te] <- mean(ytr)
      next
    }
    Xtr <- X_all[tr, keep_outer, drop = FALSE]
    Xte <- X_all[te, keep_outer, drop = FALSE]
    feat_outer <- colnames(Xtr)

    for (std_flag in c(FALSE, TRUE)) {
      sfx <- if (std_flag) "_std" else "_nostd"
      sd1 <- if (std_flag) 1 else 0
      st  <- prep_std(Xtr, Xte, std_flag)

      # L1 (lasso) and L2 (ridge) regularised linear regression
      for (alpha_g in c(1, 0)) {
        tag <- paste0(if (alpha_g == 1) "lasso" else "ridge", sfx)
        tune <- tune_glmnet_innerK(Xtr, ytr, innerK, std_flag, alpha_g,
                                   seed = seed + 1000 * i + 10 * alpha_g + sd1)
        fit <- glmnet::glmnet(x = Xtr, y = ytr, family = "gaussian", alpha = alpha_g,
                              lambda = tune$lambda, standardize = std_flag)
        pred_store[[tag]][te] <- as.numeric(predict(fit, newx = Xte))
        store_imp(tag, i, importance_glmnet(fit, tune$lambda, feat_outer), feat_outer)
      }

      # linear SVR
      tag <- paste0("svr_lin", sfx)
      tune <- tune_svr_innerK(Xtr, ytr, innerK, std_flag, svr_cost_grid, kernel = "linear",
                              seed = seed + 2000 * i + sd1)
      fit <- fit_svr(st$Xtr, ytr, "linear", tune$cost)
      pred_store[[tag]][te] <- as.numeric(predict(fit, st$Xva))
      store_imp(tag, i, importance_svr_linear(fit, st$Xtr, st$sigma), feat_outer)

      # polynomial-kernel SVR
      tag <- paste0("svr_poly", sfx)
      tune <- tune_svr_innerK(Xtr, ytr, innerK, std_flag, svr_poly_cost_grid, kernel = "polynomial",
                              degree_grid = svr_poly_degree_grid, seed = seed + 3000 * i + sd1)
      fit <- fit_svr(st$Xtr, ytr, "polynomial", tune$cost, tune$degree)
      pred_store[[tag]][te] <- as.numeric(predict(fit, st$Xva))
      store_imp(tag, i, importance_perm_rho(fit, st$Xtr, ytr, nperm = svr_poly_nperm,
                                            seed = seed + 3100 * i + sd1), feat_outer)

      # regression tree
      tag <- paste0("cart", sfx)
      tune <- tune_cart_innerK(Xtr, ytr, innerK, std_flag, cart_cp_grid, seed = seed + 5000 * i + sd1)
      fit <- fit_cart(st$Xtr, ytr, tune$cp)
      pred_store[[tag]][te] <- as.numeric(predict(fit, newdata = as.data.frame(st$Xva)))
      store_imp(tag, i, importance_cart_vi(fit, feat_outer), feat_outer)

      # gradient boosting
      tag <- paste0("gbm", sfx)
      tune <- tune_gbm_innerK(Xtr, ytr, innerK, std_flag, gbm_depth_grid, gbm_trees_grid,
                              seed = seed + 6000 * i + sd1)
      fit <- fit_gbm(st$Xtr, ytr, tune$depth, seed = seed + 7777 * i + sd1)
      pred_store[[tag]][te] <- as.numeric(predict(fit, newdata = as.data.frame(st$Xva), n.trees = tune$trees))
      store_imp(tag, i, importance_gbm_rel(fit, tune$trees, feat_outer), feat_outer)
    }
  }

  # per-accession predictions
  fold_results <- data.frame(accession = ids[tested_idx], code = y_all[tested_idx],
                             set = ifelse(tested_idx %in% train_idx, "train27", "test150"),
                             has_unseen_haplotype = has_unseen_level[tested_idx],
                             stringsAsFactors = FALSE)
  for (nm in model_names) fold_results[[paste0(nm, "_pred")]] <- round(pred_store[[nm]][tested_idx], 4)
  write.csv(fold_results, file.path(run_dir, "fold_results.csv"), row.names = FALSE)

  # metrics over the tested rows; for holdout150 also over the accessions whose
  # haplotypes were all present in the 27
  metric_row <- function(pred_vec, idx) {
    if (length(idx) < 5) return(c(rho = NA_real_, p = NA_real_, mae = NA_real_,
                                  rmse = NA_real_, mae_baseline = NA_real_))
    p <- pred_vec[idx]; y <- y_all[idx]
    c(rho = safe_cor(p, y), p = safe_cor_p(p, y), mae = mean(abs(p - y)),
      rmse = sqrt(mean((p - y)^2)), mae_baseline = mean(abs(train_median[idx] - y)))
  }
  known_idx <- if (regime == "holdout150") tested_idx[!has_unseen_level[tested_idx]] else tested_idx
  do.call(rbind, lapply(model_names, function(nm) {
    write.csv(summarise_importance(IMP[[nm]], signed = nm %in% signed_models, label_of = label_of),
              file.path(run_dir, paste0(nm, "_importance.csv")), row.names = FALSE)
    m <- metric_row(pred_store[[nm]], tested_idx); mk <- metric_row(pred_store[[nm]], known_idx)
    data.frame(model = nm, seed = seed, regime = regime, n_tested = length(tested_idx),
               spearman = m[["rho"]], spearman_p = m[["p"]], mae = m[["mae"]], rmse = m[["rmse"]],
               mae_baseline = m[["mae_baseline"]],
               n_known_haplotypes = length(known_idx), spearman_known_haplotypes = mk[["rho"]],
               stringsAsFactors = FALSE)
  }))
}

# =============================================================================
# RUN ALL REGIMES AND SEEDS
# =============================================================================
all_means <- list()
for (regime in regimes_to_run) {
  cat("\n== regime:", regime, "==\n")
  regime_dir <- file.path(out_dir, regime)
  comp_all <- dplyr::bind_rows(lapply(seeds_to_run, function(s) {
    res <- run_one_seed(seed = s, regime = regime, run_dir = file.path(regime_dir, paste0("seed_", s)))
    cat("  seed", s, "done\n")
    res
  }))
  mean_metrics <- comp_all %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      mean_spearman = mean(spearman), sd_spearman = sd(spearman),
      median_spearman_p = median(spearman_p, na.rm = TRUE),
      mean_mae = mean(mae), mean_rmse = mean(rmse), mean_mae_baseline = mean(mae_baseline),
      mean_spearman_known_haplotypes = mean(spearman_known_haplotypes, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_spearman))
  write.csv(mean_metrics, file.path(regime_dir, "mean_metrics_over_seeds.csv"), row.names = FALSE)
  print(as.data.frame(mean_metrics[, c("model", "mean_spearman", "median_spearman_p", "mean_mae",
                                       "mean_mae_baseline")]), row.names = FALSE, digits = 3)
  mean_metrics$regime <- regime
  all_means[[regime]] <- mean_metrics
}
metrics_all <- dplyr::bind_rows(all_means)
write.csv(metrics_all, file.path(out_dir, "METRICS_by_model_and_regime.csv"), row.names = FALSE)
write.csv(tidyr::pivot_wider(metrics_all[, c("model", "regime", "mean_spearman")],
                             names_from = regime, values_from = mean_spearman),
          file.path(out_dir, "METRICS_spearman_wide.csv"), row.names = FALSE)

# =============================================================================
# IMPORTANCE CONSENSUS ACROSS MODELS
# =============================================================================
# Within each retained model and seed, haplotypes are ranked by importance
# magnitude (1 = largest); haplotypes with zero importance are left unranked.
# Ranks are averaged over seeds within a model, then summarised across models
# by the median, the number of models ranking the haplotype and, from the
# signed models that ranked it, the majority direction of effect.
if (importance_regime %in% regimes_to_run) {
  mm <- all_means[[importance_regime]]
  keep <- mm$model[!is.na(mm$mean_spearman) & mm$mean_spearman > importance_threshold]
  cat("\nImportance consensus over", length(keep), "of", nrow(mm), "models with Spearman >",
      importance_threshold, "in", importance_regime, "\n")
  if (length(keep) < nrow(mm)) cat("  excluded:", paste(setdiff(mm$model, keep), collapse = ", "), "\n")
  if (!length(keep)) {
    cat("  no model cleared the threshold: no consensus written\n")
  } else {
    long <- dplyr::bind_rows(lapply(keep, function(m) {
      per_seed <- lapply(seeds_to_run, function(s) {
        d <- read.csv(file.path(out_dir, importance_regime, paste0("seed_", s),
                                paste0(m, "_importance.csv")), stringsAsFactors = FALSE)
        if (length(unique(d$magnitude)) < 2) return(NULL)   # no ranking information
        d$rank <- rank(-d$magnitude, ties.method = "min")
        d$rank[d$magnitude <= 0] <- NA_integer_
        if (!"direction" %in% names(d)) d$direction <- NA_character_
        d$direction[is.na(d$rank)] <- NA_character_   # no direction without a weight
        d[, c("haplotype", "rank", "direction")]
      })
      dd <- dplyr::bind_rows(per_seed)
      if (!nrow(dd)) return(NULL)
      dd %>%
        dplyr::group_by(haplotype) %>%
        dplyr::summarise(
          rank = if (all(is.na(rank))) NA_real_ else mean(rank, na.rm = TRUE),
          direction = if (all(is.na(direction))) NA_character_ else
            names(sort(table(direction), decreasing = TRUE))[1],
          .groups = "drop") %>%
        dplyr::mutate(model = m, model_score = mm$mean_spearman[match(m, mm$model)])
    }))
    write.csv(long, file.path(out_dir, paste0("IMPORTANCE_ranks_long_", importance_regime, ".csv")),
              row.names = FALSE)
    consensus <- long %>%
      dplyr::group_by(haplotype) %>%
      dplyr::summarise(
        n_models_using_it = sum(!is.na(rank)),
        median_rank       = median(rank, na.rm = TRUE),
        pos_votes         = sum(!is.na(rank) & !is.na(direction) & direction == dir_high),
        neg_votes         = sum(!is.na(rank) & !is.na(direction) & direction == dir_low),
        .groups = "drop") %>%
      dplyr::mutate(consensus_direction = dplyr::case_when(
        pos_votes > neg_votes ~ dir_high, neg_votes > pos_votes ~ dir_low, TRUE ~ "mixed/none")) %>%
      dplyr::arrange(median_rank, dplyr::desc(n_models_using_it))
    write.csv(consensus, file.path(out_dir, paste0("IMPORTANCE_consensus_", importance_regime, ".csv")),
              row.names = FALSE)
    cat("\nTop 15 haplotypes by median rank across the", length(unique(long$model)), "retained models:\n")
    print(as.data.frame(head(consensus, 15)), row.names = FALSE, digits = 3)
  }
}

# =============================================================================
# ASSOCIATION OF EVERY HAPLOTYPE WITH THE SCORE, ALL 177 ACCESSIONS
# =============================================================================
very <- y_all >= very_compact_from
assoc <- dplyr::bind_rows(lapply(haplo_cols, function(m) {
  vals <- as.character(df0[[m]])
  dplyr::bind_rows(lapply(sort(unique(vals)), function(hp) {
    has <- vals == hp; n <- sum(has)
    tested <- n >= assoc_min_carriers && n <= length(y_all) - assoc_min_carriers
    ft <- if (tested) stats::fisher.test(matrix(c(sum(has & very), n - sum(has & very),
                                                  sum(!has & very), sum(!has) - sum(!has & very)), 2)) else NULL
    data.frame(
      locus = m, haplotype = hp, label = paste0(m, ":", hp), n_carriers = n,
      mean_code_carriers = mean(y_all[has]), mean_code_others = mean(y_all[!has]),
      mean_code_diff = mean(y_all[has]) - mean(y_all[!has]),
      p_mannwhitney = if (tested)
        suppressWarnings(stats::wilcox.test(y_all[has], y_all[!has], exact = FALSE)$p.value) else NA_real_,
      n_very_compact_carriers = sum(has & very),
      share_very_compact_carriers = sum(has & very) / n,
      odds_ratio_very_compact = if (tested) unname(ft$estimate) else NA_real_,
      p_fisher_very_compact = if (tested) ft$p.value else NA_real_,
      stringsAsFactors = FALSE)
  }))
}))
ok <- !is.na(assoc$p_mannwhitney)
assoc$q_mannwhitney <- NA_real_; assoc$q_fisher_very_compact <- NA_real_
assoc$q_mannwhitney[ok] <- stats::p.adjust(assoc$p_mannwhitney[ok], method = "BH")
assoc$q_fisher_very_compact[ok] <- stats::p.adjust(assoc$p_fisher_very_compact[ok], method = "BH")
assoc$direction <- ifelse(assoc$mean_code_diff > 0, dir_high, dir_low)
assoc$significant <- !is.na(assoc$q_mannwhitney) & assoc$q_mannwhitney < assoc_alpha
assoc <- assoc[order(assoc$q_mannwhitney, -abs(assoc$mean_code_diff), na.last = TRUE), ]
write.csv(assoc, file.path(out_dir, "HAPLOTYPE_association.csv"), row.names = FALSE)

cat(sprintf("\nAssociation with the score on all %d accessions (Mann-Whitney, BH q < %.2f): %d of %d tested haplotypes\n",
            length(y_all), assoc_alpha, sum(assoc$significant), sum(ok)))
print(head(assoc[ok, c("label", "n_carriers", "mean_code_carriers", "mean_code_others",
                       "q_mannwhitney", "direction")], 10), row.names = FALSE, digits = 3)
cat(sprintf("Very compact panicles (code >= %d, n = %d), haplotypes at q < %.2f: %s\n",
            very_compact_from, sum(very), assoc_alpha,
            paste(assoc$label[!is.na(assoc$q_fisher_very_compact) &
                              assoc$q_fisher_very_compact < assoc_alpha], collapse = ", ")))

cat("\nOutputs written under:", normalizePath(out_dir), "\n")
