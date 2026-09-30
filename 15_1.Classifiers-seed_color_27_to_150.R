###############################################################################
# Haplotype-based classification of teff seed colour (White vs Brown)
# ---------------------------------------------------------------------------
# Fourteen classifier variants are trained on the 27 pangenome accessions and
# evaluated (i) internally, on the 27, by resubstitution, leave-one-out
# cross-validation and stratified 5-fold cross-validation, and (ii) externally,
# on the 150 EIAR accessions that are never used for training or tuning.
# Haplotype importance is then summarised as a rank consensus across the
# variants that classify the 150 well.
#
# INPUT
#   new_data.txt   tab-separated, one row per accession: Name, phenotype
#                  (White / Brown) and one "haplotype-<locus>" column per locus,
#                  with haplotypes called jointly on all 177 accessions.
#                  Loci 4B and 4C are the two intervals on chromosome 4B
#                  (4B1 and 4B2 in the manuscript); 6A and 6B are the two
#                  intervals on chromosome 6A (6A1 and 6A2).
#
# MODELS (each fitted with and without predictor standardisation)
#   L1- and L2-regularised logistic regression (glmnet), linear and
#   polynomial-kernel support vector machines (e1071), naive Bayes (e1071),
#   classification trees (rpart) and gradient boosting (gbm). Haplotypes enter
#   as binary indicators (one per haplotype at each locus). Hyperparameters are
#   tuned by stratified inner cross-validation on the training rows of each
#   fold, with balanced accuracy as the criterion; standardisation uses
#   training-fold statistics only; indicators constant in the training rows of
#   a fold are removed before fitting.
#
# EVALUATION REGIMES
#   resub27     fitted and scored on the 27 (in-sample fit)
#   loocv27     leave-one-out cross-validation on the 27
#   cv5_27      stratified 5-fold cross-validation on the 27
#   holdout150  fitted on all 27, scored on the 150 EIAR accessions
# Each regime is repeated for every seed in seeds_to_run and averaged.
#
# UNSEEN HAPLOTYPES
#   A haplotype present among the 150 but absent from the 27 has no indicator
#   in a model fitted on the 27, so accessions carrying it contribute zeros at
#   that locus. They are kept and flagged, and holdout150 metrics are reported
#   both on all 150 and on the accessions carrying only training-set
#   haplotypes.
#
# OUTPUTS (under out_dir)
#   accessions_split.csv                  set (train27 / test150), phenotype,
#                                         unseen-haplotype flag, per accession
#   unseen_haplotypes.csv                 per locus: haplotypes present in the
#                                         27, and number of test accessions
#                                         carrying a haplotype absent from the 27
#   <regime>/seed_<s>/fold_results.csv    per-accession predictions
#   <regime>/seed_<s>/<model>_importance.csv
#                                         per-haplotype importance of one model
#   <regime>/mean_metrics_over_seeds.csv  metrics averaged over seeds
#   METRICS_by_model_and_regime.csv       all regimes, means and sd over seeds
#   METRICS_accuracy_wide.csv             model x regime, mean accuracy
#   METRICS_balanced_accuracy_wide.csv    model x regime, mean balanced accuracy
#   IMPORTANCE_ranks_long_holdout150.csv  rank of every haplotype in every
#                                         retained model
#   IMPORTANCE_consensus_holdout150.csv   median rank across models, number of
#                                         models ranking the haplotype,
#                                         direction from the signed models
#   MISCLASSIFIED_by_all_models_holdout150.csv
#                                         test accessions that no model
#                                         classifies correctly in any seed
#
# Runtime: a few minutes per seed and regime (training sets are 27 rows).
###############################################################################

# =============================================================================
# CONFIGURATION
# =============================================================================
data_file <- "new_data.txt"
out_dir   <- "outputs_classifiers"

# The accession recorded as "DZ" in new_data.txt is DZ_01_354 in the manuscript.
accession_rename <- c(DZ = "DZ_01_354")

# The 27 pangenome accessions (training set). All other rows are the 150 EIAR
# accessions used for external validation.
train_names_27 <- c(
  "Boni", "DZ_01_354", "Dabbi", "Dtt2-02", "Quncho", "T-33", "T116", "T132",
  "T177", "T206", "T224", "T283", "T288", "T297", "T304", "T330",
  "T336", "T345", "T365", "T366", "T379", "T404", "T412", "T87",
  "T99", "addisie", "karadebi"
)

seeds_to_run   <- 1:3
regimes_to_run <- c("resub27", "loocv27", "cv5_27", "holdout150")
cv_folds       <- 5      # outer folds of the cv5_27 regime
innerK_base    <- 5      # inner folds for hyperparameter tuning (reduced if a
                         # class has fewer training accessions than this)

# Hyperparameter grids
svm_cost_grid        <- c(0.1, 1, 10, 100)
svm_poly_cost_grid   <- c(0.1, 1, 10)
svm_poly_degree_grid <- c(2, 3)
svm_poly_coef0       <- 0
nb_prior_grid        <- c("empirical", "uniform")
nb_laplace_grid      <- c(0, 1)
cart_cp_grid         <- c(0.0, 0.001, 0.01, 0.05)
cart_maxdepth        <- 30
cart_minsplit        <- 2
gbm_depth_grid       <- c(1, 2, 3)
gbm_trees_grid       <- c(50, 100, 200, 500)
gbm_shrinkage        <- 0.01
gbm_minobs           <- 2
gbm_bagfrac          <- 0.8
svm_poly_nperm       <- 10   # permutations for the polynomial-SVM importance

# Importance consensus: models are retained when their mean balanced accuracy
# in importance_regime exceeds importance_threshold.
importance_regime    <- "holdout150"
importance_metric    <- "balanced_accuracy"   # or "accuracy"
importance_threshold <- 0.75

model_names <- c(
  "lasso_nostd", "lasso_std", "l2log_nostd", "l2log_std",
  "svm_lin_nostd", "svm_lin_std", "svm_poly_nostd", "svm_poly_std",
  "nb_nostd", "nb_std", "cart_nostd", "cart_std", "gbm_nostd", "gbm_std"
)
signed_models <- c("lasso_nostd", "lasso_std", "l2log_nostd", "l2log_std",
                   "svm_lin_nostd", "svm_lin_std")

pkgs <- c("dplyr", "tidyr", "glmnet", "e1071", "rpart", "gbm")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# METRICS AND HELPERS
# =============================================================================
clf_metrics <- function(y_true, y_pred, positive_class = "Brown") {
  y_true <- droplevels(factor(y_true))
  y_pred <- factor(y_pred, levels = levels(y_true))
  if (nlevels(y_true) < 2) {
    return(list(acc = mean(y_true == y_pred), bal_acc = NA_real_,
                sensitivity = NA_real_, specificity = NA_real_))
  }
  neg_class <- setdiff(levels(y_true), positive_class)[1]
  TP <- sum(y_true == positive_class & y_pred == positive_class)
  TN <- sum(y_true == neg_class      & y_pred == neg_class)
  FP <- sum(y_true == neg_class      & y_pred == positive_class)
  FN <- sum(y_true == positive_class & y_pred == neg_class)
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  list(acc = (TP + TN) / (TP + TN + FP + FN),
       bal_acc = mean(c(sens, spec), na.rm = TRUE),
       sensitivity = sens, specificity = spec)
}

bal_acc_safe <- function(y_true, y_pred) {
  m <- clf_metrics(y_true, y_pred, positive_class = "Brown")$bal_acc
  if (is.na(m)) 0 else m
}

# Standardise with training statistics only.
standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sigma <- apply(Xtr, 2, sd)
  sigma[is.na(sigma) | sigma == 0] <- 1
  list(Xtr = sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/"),
       Xte = sweep(sweep(Xte, 2, mu, "-"), 2, sigma, "/"),
       mu = mu, sigma = sigma)
}

prep_std <- function(Xtr, Xva, standardize_flag) {
  if (!standardize_flag) return(list(Xtr = Xtr, Xva = Xva, sigma = rep(1, ncol(Xtr))))
  st <- standardize_train_test(Xtr, Xva)
  list(Xtr = st$Xtr, Xva = st$Xte, sigma = st$sigma)
}

make_stratified_folds <- function(y, K = 5, seed = 1) {
  set.seed(seed)
  y <- droplevels(factor(y))
  fold_id <- rep(NA_integer_, length(y))
  for (cl in levels(y)) {
    idx <- which(y == cl)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(1:K, length.out = length(idx))
  }
  fold_id
}

# =============================================================================
# INNER-CV TUNERS (criterion: balanced accuracy)
# =============================================================================
tune_glmnet_innerK <- function(Xtr, ytr_bin, ytr_fac, K, standardize_flag,
                               alpha_glmnet, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  fit0 <- glmnet::glmnet(x = Xtr, y = ytr_bin, family = "binomial",
                         alpha = alpha_glmnet, standardize = standardize_flag,
                         intercept = FALSE)
  lambda_grid <- fit0$lambda
  fold_scores <- matrix(NA_real_, nrow = K, ncol = length(lambda_grid))
  for (k in 1:K) {
    va <- which(folds == k); tr <- which(folds != k)
    if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_bin[tr])) < 2) next
    fit <- glmnet::glmnet(x = Xtr[tr, , drop = FALSE], y = ytr_bin[tr],
                          family = "binomial", alpha = alpha_glmnet,
                          lambda = lambda_grid, standardize = standardize_flag,
                          intercept = FALSE)
    pmat <- predict(fit, newx = Xtr[va, , drop = FALSE], type = "response")
    idx_used <- match(fit$lambda, lambda_grid)
    for (j in seq_along(fit$lambda)) {
      pred <- ifelse(as.numeric(pmat[, j]) >= 0.5, "Brown", "White")
      fold_scores[k, idx_used[j]] <- bal_acc_safe(ytr_fac[va], pred)
    }
  }
  # one-standard-error rule: the most regularised lambda within one SE of the best
  mean_sc <- colMeans(fold_scores, na.rm = TRUE)
  sd_sc   <- apply(fold_scores, 2, sd, na.rm = TRUE)
  se_sc   <- sd_sc / sqrt(pmax(colSums(!is.na(fold_scores)), 1))
  j_best  <- which.max(mean_sc)
  cand    <- which(mean_sc >= max(mean_sc, na.rm = TRUE) - se_sc[j_best])
  if (length(cand) == 0) cand <- j_best
  list(lambda = lambda_grid[max(cand)])
}

tune_svm_linear_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                                   cost_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  scores <- setNames(rep(NA_real_, length(cost_grid)), as.character(cost_grid))
  for (cst in cost_grid) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- e1071::svm(x = st$Xtr, y = droplevels(ytr_fac[tr]),
                        kernel = "linear", cost = cst, scale = FALSE)
      fold_sc[k] <- bal_acc_safe(ytr_fac[va], predict(fit, st$Xva))
    }
    scores[as.character(cst)] <- mean(fold_sc, na.rm = TRUE)
  }
  list(cost = as.numeric(names(which.max(scores))))
}

tune_svm_poly_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                                 cost_grid, degree_grid, coef0, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(cost = cost_grid, degree = degree_grid, stringsAsFactors = FALSE)
  scores <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- e1071::svm(x = st$Xtr, y = droplevels(ytr_fac[tr]),
                        kernel = "polynomial", cost = grid$cost[g],
                        degree = grid$degree[g], coef0 = coef0, scale = FALSE)
      fold_sc[k] <- bal_acc_safe(ytr_fac[va], predict(fit, st$Xva))
    }
    scores[g] <- mean(fold_sc, na.rm = TRUE)
  }
  best <- grid[which.max(scores), , drop = FALSE]
  list(cost = best$cost, degree = best$degree, coef0 = coef0)
}

fit_nb <- function(X, y, prior_type, laplace) {
  if (prior_type == "uniform") {
    pr0 <- setNames(rep(1 / nlevels(y), nlevels(y)), levels(y))
    e1071::naiveBayes(x = X, y = y, prior = pr0, laplace = laplace)
  } else {
    e1071::naiveBayes(x = X, y = y, laplace = laplace)
  }
}

tune_nb_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                           prior_grid, laplace_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(prior_type = prior_grid, laplace = laplace_grid,
                      stringsAsFactors = FALSE)
  scores <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_nb(st$Xtr, droplevels(ytr_fac[tr]), grid$prior_type[g], grid$laplace[g])
      fold_sc[k] <- bal_acc_safe(ytr_fac[va], predict(fit, st$Xva))
    }
    scores[g] <- mean(fold_sc, na.rm = TRUE)
  }
  best <- which.max(scores)
  list(prior_type = grid$prior_type[best], laplace = grid$laplace[best])
}

fit_cart <- function(X, y, cp) {
  df_tr <- data.frame(y = y, as.data.frame(X))
  ctrl <- rpart::rpart.control(cp = cp, maxdepth = cart_maxdepth,
                               minsplit = cart_minsplit, xval = 0)
  rpart::rpart(y ~ ., data = df_tr, method = "class", control = ctrl)
}

tune_cart_innerK <- function(Xtr, ytr_fac, K, standardize_flag, cp_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  scores <- setNames(rep(NA_real_, length(cp_grid)), as.character(cp_grid))
  for (cpv in cp_grid) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_cart(st$Xtr, droplevels(ytr_fac[tr]), cpv)
      pr  <- predict(fit, newdata = as.data.frame(st$Xva), type = "class")
      fold_sc[k] <- bal_acc_safe(ytr_fac[va], pr)
    }
    scores[as.character(cpv)] <- mean(fold_sc, na.rm = TRUE)
  }
  list(cp = as.numeric(names(which.max(scores))))
}

fit_gbm <- function(X, y_bin, depth, seed) {
  df_tr <- data.frame(y = y_bin, as.data.frame(X))
  set.seed(seed)
  gbm::gbm(y ~ ., data = df_tr, distribution = "bernoulli",
           n.trees = max(gbm_trees_grid), interaction.depth = depth,
           shrinkage = gbm_shrinkage, n.minobsinnode = gbm_minobs,
           bag.fraction = gbm_bagfrac, train.fraction = 1.0, cv.folds = 0,
           verbose = FALSE)
}

tune_gbm_innerK <- function(Xtr, ytr_bin, ytr_fac, K, standardize_flag,
                            depth_grid, trees_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(depth = depth_grid, trees = trees_grid, stringsAsFactors = FALSE)
  scores <- rep(NA_real_, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    for (k in 1:K) {
      va <- which(folds == k); tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      st <- prep_std(Xtr[tr, , drop = FALSE], Xtr[va, , drop = FALSE], standardize_flag)
      fit <- fit_gbm(st$Xtr, ytr_bin[tr], grid$depth[g],
                     seed = seed + 10000 * g + 100 * k + if (standardize_flag) 1 else 0)
      p_va <- predict(fit, newdata = as.data.frame(st$Xva),
                      n.trees = grid$trees[g], type = "response")
      fold_sc[k] <- bal_acc_safe(ytr_fac[va], ifelse(p_va >= 0.5, "Brown", "White"))
    }
    scores[g] <- mean(fold_sc, na.rm = TRUE)
  }
  best <- which.max(scores)
  list(depth = grid$depth[best], trees = grid$trees[best])
}

# =============================================================================
# IMPORTANCE EXTRACTORS
# =============================================================================
# Signed (positive = towards Brown): coefficients of the logistic models and
# weights of the linear SVM. Unsigned: permutation drop in balanced accuracy
# (polynomial SVM), standardised difference of indicator means between classes
# (naive Bayes), variable importance (CART) and relative influence (GBM).
importance_glmnet <- function(fit, lambda, feat_names) {
  b <- as.matrix(coef(fit, s = lambda))
  b <- b[setdiff(rownames(b), "(Intercept)"), , drop = FALSE]
  v <- setNames(rep(0, length(feat_names)), feat_names)
  common <- intersect(rownames(b), feat_names)
  v[common] <- b[common, 1]
  v
}

importance_svm_linear <- function(svm_fit, Xtr_fit, ytr, sigma, positive_class = "Brown") {
  w_scaled <- drop(t(svm_fit$coefs) %*% svm_fit$SV)
  names(w_scaled) <- colnames(Xtr_fit)
  # orient the weight vector so that positive values point towards Brown
  dv_tr <- as.numeric(attr(predict(svm_fit, Xtr_fit, decision.values = TRUE), "decision.values"))
  sgn <- ifelse(mean(dv_tr[ytr == positive_class]) - mean(dv_tr[ytr != positive_class]) >= 0, 1, -1)
  (w_scaled / sigma) * sgn        # back to the unstandardised scale
}

importance_svm_poly_perm <- function(fit, Xtr_fit, ytr, nperm = 10, seed = 1) {
  set.seed(seed)
  base_bal <- bal_acc_safe(ytr, predict(fit, Xtr_fit))
  imp <- setNames(rep(0, ncol(Xtr_fit)), colnames(Xtr_fit))
  for (j in seq_len(ncol(Xtr_fit))) {
    drops <- numeric(nperm)
    for (b in 1:nperm) {
      Xp <- Xtr_fit
      Xp[, j] <- sample(Xp[, j], replace = FALSE)
      drops[b] <- base_bal - bal_acc_safe(ytr, predict(fit, Xp))
    }
    imp[j] <- mean(pmax(drops, 0))
  }
  imp
}

importance_nb_effect <- function(Xtr_fit, ytr, positive_class = "Brown") {
  Xpos <- Xtr_fit[ytr == positive_class, , drop = FALSE]
  Xneg <- Xtr_fit[ytr != positive_class, , drop = FALSE]
  sd_pos <- apply(Xpos, 2, sd); sd_pos[is.na(sd_pos)] <- 0
  sd_neg <- apply(Xneg, 2, sd); sd_neg[is.na(sd_neg)] <- 0
  abs(colMeans(Xpos) - colMeans(Xneg)) / (sqrt((sd_pos^2 + sd_neg^2) / 2) + 1e-8)
}

importance_cart_vi <- function(fit, feat_names) {
  v <- setNames(rep(0, length(feat_names)), feat_names)
  vi <- fit$variable.importance
  if (!is.null(vi) && length(vi) > 0) {
    common <- intersect(names(vi), feat_names)
    v[common] <- as.numeric(vi[common])
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

importance_gbm_rel <- function(fit, n.trees, feat_names) {
  v <- setNames(rep(0, length(feat_names)), feat_names)
  rel <- suppressWarnings(summary(fit, plotit = FALSE, n.trees = n.trees))
  if (!is.null(rel) && nrow(rel) > 0) {
    common <- intersect(as.character(rel$var), feat_names)
    v[common] <- rel$rel.inf[match(common, as.character(rel$var))]
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

# Fold-wise importances of one model -> one row per haplotype. `magnitude` is
# the absolute mean across folds; `direction` (signed models) is the sign of
# that mean.
summarise_importance <- function(M, signed, label_of) {
  mean_i <- colMeans(M, na.rm = TRUE)
  out <- data.frame(feature = colnames(M),
                    haplotype = unname(label_of[colnames(M)]),
                    mean = as.numeric(mean_i),
                    sd = as.numeric(apply(M, 2, sd, na.rm = TRUE)),
                    magnitude = abs(as.numeric(mean_i)),
                    stringsAsFactors = FALSE)
  if (signed) out$direction <- ifelse(mean_i > 0, "Brown", "White")
  out[order(-out$magnitude), ]
}

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

hap_raw    <- setdiff(names(raw), c("Name", "phenotype"))
haplo_cols <- paste("Chr", sub("^haplotype[-_]?", "", hap_raw))     # "haplotype-4B" -> "Chr 4B"
df0 <- data.frame(phenotype = factor(raw$phenotype, levels = c("White", "Brown")))
for (j in seq_along(hap_raw)) df0[[haplo_cols[j]]] <- factor(as.character(raw[[hap_raw[j]]]))
rownames(df0) <- raw$Name

train_idx <- match(train_names_27, raw$Name)
if (anyNA(train_idx)) stop("Training accessions not found in ", data_file, ": ",
                           paste(train_names_27[is.na(train_idx)], collapse = ", "))
test_idx <- setdiff(seq_len(nrow(df0)), train_idx)

# One binary indicator per haplotype at each locus, over all accessions.
# Column names are syntactically valid R names (needed by rpart and gbm);
# `feat_label` carries the display form locus:haplotype, e.g. "4B:3".
X_all <- do.call(cbind, lapply(haplo_cols, function(m) {
  lv  <- levels(df0[[m]])
  mat <- matrix(as.numeric(outer(as.character(df0[[m]]), lv, "==")), ncol = length(lv))
  colnames(mat) <- make.names(paste0(m, "_", lv))
  mat
}))
rownames(X_all) <- rownames(df0)
feat_names   <- colnames(X_all)
feat_label   <- unlist(lapply(haplo_cols, function(m)
  paste0(sub("^Chr ", "", m), ":", levels(df0[[m]]))))
label_of     <- setNames(feat_label, feat_names)
locus_of_col <- rep(haplo_cols, vapply(haplo_cols, function(m) nlevels(df0[[m]]), integer(1)))
y_fac <- df0$phenotype
y_bin <- ifelse(y_fac == "Brown", 1, 0)
ids   <- rownames(df0)

# Haplotypes absent from the 27: a test accession has an unseen haplotype at a
# locus when it carries none of that locus's haplotypes present in the 27.
present_in_27 <- colSums(X_all[train_idx, , drop = FALSE]) > 0
unseen_mat <- vapply(haplo_cols, function(m) {
  cols <- which(locus_of_col == m & present_in_27)
  rowSums(X_all[test_idx, cols, drop = FALSE]) == 0
}, logical(length(test_idx)))
has_unseen_level <- rep(NA, nrow(df0))
has_unseen_level[test_idx] <- rowSums(unseen_mat) > 0

cat("Training accessions:", length(train_idx), "| validation accessions:", length(test_idx),
    "| loci:", length(haplo_cols), "| haplotype indicators:", ncol(X_all), "\n")
cat("Validation accessions carrying at least one haplotype absent from the 27:",
    sum(has_unseen_level[test_idx]), "\n")

write.csv(data.frame(accession = ids,
                     set = ifelse(seq_along(ids) %in% train_idx, "train27", "test150"),
                     phenotype = as.character(y_fac),
                     has_unseen_haplotype = has_unseen_level),
          file.path(out_dir, "accessions_split.csv"), row.names = FALSE)
write.csv(data.frame(locus = haplo_cols,
                     haplotypes_in_27 = vapply(haplo_cols, function(m)
                       paste(levels(droplevels(df0[[m]][train_idx])), collapse = "/"), ""),
                     n_test_with_unseen_haplotype = as.integer(colSums(unseen_mat))),
          file.path(out_dir, "unseen_haplotypes.csv"), row.names = FALSE)

# =============================================================================
# FOLD PLANS: every regime is a list of (train, test) row-index pairs
# =============================================================================
make_fold_plan <- function(regime, seed) {
  tr <- train_idx
  switch(regime,
    resub27    = list(list(train = tr, test = tr)),
    loocv27    = lapply(seq_along(tr), function(j) list(train = tr[-j], test = tr[j])),
    cv5_27     = {
      fid <- make_stratified_folds(y_fac[tr], K = cv_folds, seed = seed)
      lapply(sort(unique(fid)), function(k) list(train = tr[fid != k], test = tr[fid == k]))
    },
    holdout150 = list(list(train = tr, test = test_idx)),
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

  pred_store <- setNames(lapply(model_names, function(x) rep(NA_character_, n)), model_names)
  IMP <- setNames(lapply(model_names, function(x)
    matrix(0, nrow = n_folds, ncol = length(feat_names),
           dimnames = list(paste0("fold", seq_len(n_folds)), feat_names))), model_names)
  store_imp <- function(tag, i, v_keep, feat_outer) {
    v_full <- setNames(rep(0, length(feat_names)), feat_names)
    v_full[feat_outer] <- v_keep
    IMP[[tag]][i, ] <<- v_full
  }

  for (i in seq_len(n_folds)) {
    te <- fold_plan[[i]]$test
    tr <- fold_plan[[i]]$train
    ytr_fac <- y_fac[tr]; ytr_bin <- y_bin[tr]
    innerK  <- max(2, min(innerK_base, min(table(droplevels(ytr_fac)))))

    # indicators constant in the training rows carry no information for this fold
    sd_outer <- apply(X_all[tr, , drop = FALSE], 2, sd); sd_outer[is.na(sd_outer)] <- 0
    keep_outer <- sd_outer > 0
    if (!any(keep_outer)) {
      maj <- names(which.max(table(ytr_fac)))
      for (nm in model_names) pred_store[[nm]][te] <- maj
      next
    }
    Xtr <- X_all[tr, keep_outer, drop = FALSE]
    Xte <- X_all[te, keep_outer, drop = FALSE]
    feat_outer <- colnames(Xtr)
    ytr_drop <- droplevels(ytr_fac)

    for (std_flag in c(FALSE, TRUE)) {
      sfx <- if (std_flag) "_std" else "_nostd"
      sd1 <- if (std_flag) 1 else 0
      st  <- prep_std(Xtr, Xte, std_flag)

      # L1 (lasso) and L2 (ridge) regularised logistic regression
      for (alpha_glmnet in c(1, 0)) {
        tag <- paste0(if (alpha_glmnet == 1) "lasso" else "l2log", sfx)
        tune <- tune_glmnet_innerK(Xtr, ytr_bin, ytr_fac, K = innerK,
                                   standardize_flag = std_flag, alpha_glmnet = alpha_glmnet,
                                   seed = seed + (if (alpha_glmnet == 1) 1000 else 1500) * i + sd1)
        fit <- glmnet::glmnet(x = Xtr, y = ytr_bin, family = "binomial", alpha = alpha_glmnet,
                              lambda = tune$lambda, standardize = std_flag, intercept = FALSE)
        p_te <- as.numeric(predict(fit, newx = Xte, type = "response"))
        pred_store[[tag]][te] <- ifelse(p_te >= 0.5, "Brown", "White")
        store_imp(tag, i, importance_glmnet(fit, tune$lambda, feat_outer), feat_outer)
      }

      # linear SVM
      tag <- paste0("svm_lin", sfx)
      tune <- tune_svm_linear_innerK(Xtr, ytr_fac, K = innerK, standardize_flag = std_flag,
                                     cost_grid = svm_cost_grid, seed = seed + 2000 * i + sd1)
      fit <- e1071::svm(x = st$Xtr, y = ytr_drop, kernel = "linear", cost = tune$cost,
                        scale = FALSE, decision.values = TRUE)
      pred_store[[tag]][te] <- as.character(predict(fit, st$Xva, decision.values = TRUE))
      store_imp(tag, i, importance_svm_linear(fit, st$Xtr, ytr_drop, sigma = st$sigma), feat_outer)

      # polynomial-kernel SVM
      tag <- paste0("svm_poly", sfx)
      tune <- tune_svm_poly_innerK(Xtr, ytr_fac, K = innerK, standardize_flag = std_flag,
                                   cost_grid = svm_poly_cost_grid, degree_grid = svm_poly_degree_grid,
                                   coef0 = svm_poly_coef0, seed = seed + 3000 * i + sd1)
      fit <- e1071::svm(x = st$Xtr, y = ytr_drop, kernel = "polynomial", cost = tune$cost,
                        degree = tune$degree, coef0 = tune$coef0, scale = FALSE)
      pred_store[[tag]][te] <- as.character(predict(fit, st$Xva))
      store_imp(tag, i, importance_svm_poly_perm(fit, st$Xtr, ytr_drop, nperm = svm_poly_nperm,
                                                 seed = seed + 3100 * i + sd1), feat_outer)

      # naive Bayes
      tag <- paste0("nb", sfx)
      tune <- tune_nb_innerK(Xtr, ytr_fac, K = innerK, standardize_flag = std_flag,
                             prior_grid = nb_prior_grid, laplace_grid = nb_laplace_grid,
                             seed = seed + 4000 * i + sd1)
      fit <- fit_nb(st$Xtr, ytr_drop, tune$prior_type, tune$laplace)
      pred_store[[tag]][te] <- as.character(predict(fit, st$Xva))
      store_imp(tag, i, importance_nb_effect(st$Xtr, ytr_drop), feat_outer)

      # CART
      tag <- paste0("cart", sfx)
      tune <- tune_cart_innerK(Xtr, ytr_fac, K = innerK, standardize_flag = std_flag,
                               cp_grid = cart_cp_grid, seed = seed + 5000 * i + sd1)
      fit <- fit_cart(st$Xtr, ytr_drop, tune$cp)
      pred_store[[tag]][te] <- as.character(predict(fit, newdata = as.data.frame(st$Xva), type = "class"))
      store_imp(tag, i, importance_cart_vi(fit, feat_outer), feat_outer)

      # gradient boosting
      tag <- paste0("gbm", sfx)
      tune <- tune_gbm_innerK(Xtr, ytr_bin, ytr_fac, K = innerK, standardize_flag = std_flag,
                              depth_grid = gbm_depth_grid, trees_grid = gbm_trees_grid,
                              seed = seed + 6000 * i + sd1)
      fit <- fit_gbm(st$Xtr, ytr_bin, tune$depth, seed = seed + 7777 * i + sd1)
      p_te <- predict(fit, newdata = as.data.frame(st$Xva), n.trees = tune$trees, type = "response")
      pred_store[[tag]][te] <- ifelse(p_te >= 0.5, "Brown", "White")
      store_imp(tag, i, importance_gbm_rel(fit, n.trees = tune$trees, feat_names = feat_outer), feat_outer)
    }
  }

  # per-accession predictions
  truth <- as.character(y_fac)
  fold_results <- data.frame(accession = ids[tested_idx], truth = truth[tested_idx],
                             set = ifelse(tested_idx %in% train_idx, "train27", "test150"),
                             has_unseen_haplotype = has_unseen_level[tested_idx],
                             stringsAsFactors = FALSE)
  for (nm in model_names) {
    fold_results[[paste0(nm, "_pred")]]    <- pred_store[[nm]][tested_idx]
    fold_results[[paste0(nm, "_correct")]] <- as.integer(pred_store[[nm]][tested_idx] == truth[tested_idx])
  }
  write.csv(fold_results, file.path(run_dir, "fold_results.csv"), row.names = FALSE)

  # metrics over the tested rows; for holdout150 also over the accessions whose
  # haplotypes were all present in the 27
  metric_row <- function(pred_vec, idx) {
    if (length(idx) == 0 || nlevels(droplevels(y_fac[idx])) < 2)
      return(c(acc = NA_real_, bal_acc = NA_real_, sens = NA_real_, spec = NA_real_))
    met <- clf_metrics(y_fac[idx], factor(pred_vec[idx], levels = levels(y_fac)))
    c(acc = met$acc, bal_acc = met$bal_acc, sens = met$sensitivity, spec = met$specificity)
  }
  known_idx <- if (regime == "holdout150") tested_idx[!has_unseen_level[tested_idx]] else tested_idx
  comparison <- do.call(rbind, lapply(model_names, function(nm) {
    m  <- metric_row(pred_store[[nm]], tested_idx)
    mk <- metric_row(pred_store[[nm]], known_idx)
    data.frame(model = nm, seed = seed, regime = regime, n_tested = length(tested_idx),
               accuracy = m[["acc"]], balanced_accuracy = m[["bal_acc"]],
               sensitivity_brown = m[["sens"]], specificity_white = m[["spec"]],
               n_known_haplotypes = length(known_idx),
               accuracy_known_haplotypes = mk[["acc"]],
               balanced_accuracy_known_haplotypes = mk[["bal_acc"]],
               stringsAsFactors = FALSE)
  }))

  for (nm in model_names) {
    write.csv(summarise_importance(IMP[[nm]], signed = nm %in% signed_models, label_of = label_of),
              file.path(run_dir, paste0(nm, "_importance.csv")), row.names = FALSE)
  }
  comparison
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
      mean_accuracy = mean(accuracy, na.rm = TRUE), sd_accuracy = sd(accuracy, na.rm = TRUE),
      mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
      sd_balanced_accuracy = sd(balanced_accuracy, na.rm = TRUE),
      mean_sensitivity_brown = mean(sensitivity_brown, na.rm = TRUE),
      mean_specificity_white = mean(specificity_white, na.rm = TRUE),
      mean_accuracy_known_haplotypes = mean(accuracy_known_haplotypes, na.rm = TRUE),
      mean_balanced_accuracy_known_haplotypes = mean(balanced_accuracy_known_haplotypes, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_accuracy))
  write.csv(mean_metrics, file.path(regime_dir, "mean_metrics_over_seeds.csv"), row.names = FALSE)
  print(as.data.frame(mean_metrics[, c("model", "mean_accuracy", "mean_balanced_accuracy",
                                       "mean_sensitivity_brown", "mean_specificity_white")]),
        row.names = FALSE, digits = 3)
  mean_metrics$regime <- regime
  all_means[[regime]] <- mean_metrics
}
metrics_all <- dplyr::bind_rows(all_means)
write.csv(metrics_all, file.path(out_dir, "METRICS_by_model_and_regime.csv"), row.names = FALSE)
write.csv(tidyr::pivot_wider(metrics_all[, c("model", "regime", "mean_accuracy")],
                             names_from = regime, values_from = mean_accuracy),
          file.path(out_dir, "METRICS_accuracy_wide.csv"), row.names = FALSE)
write.csv(tidyr::pivot_wider(metrics_all[, c("model", "regime", "mean_balanced_accuracy")],
                             names_from = regime, values_from = mean_balanced_accuracy),
          file.path(out_dir, "METRICS_balanced_accuracy_wide.csv"), row.names = FALSE)

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
  mm$score <- mm[[paste0("mean_", importance_metric)]]
  keep <- mm$model[!is.na(mm$score) & mm$score > importance_threshold]
  cat("\nImportance consensus over", length(keep), "of", nrow(mm), "models with",
      importance_metric, ">", importance_threshold, "in", importance_regime, "\n")
  if (length(keep) < nrow(mm))
    cat("  excluded:", paste(setdiff(mm$model, keep), collapse = ", "), "\n")

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
      dplyr::mutate(model = m, model_score = mm$score[match(m, mm$model)])
  }))
  write.csv(long, file.path(out_dir, paste0("IMPORTANCE_ranks_long_", importance_regime, ".csv")),
            row.names = FALSE)

  consensus <- long %>%
    dplyr::group_by(haplotype) %>%
    dplyr::summarise(
      n_models_using_it = sum(!is.na(rank)),
      median_rank       = median(rank, na.rm = TRUE),
      brown_votes       = sum(!is.na(rank) & !is.na(direction) & direction == "Brown"),
      white_votes       = sum(!is.na(rank) & !is.na(direction) & direction == "White"),
      .groups = "drop") %>%
    dplyr::mutate(consensus_direction = dplyr::case_when(
      brown_votes > white_votes ~ "Brown",
      white_votes > brown_votes ~ "White",
      TRUE ~ "mixed/none")) %>%
    dplyr::arrange(median_rank, dplyr::desc(n_models_using_it))
  write.csv(consensus, file.path(out_dir, paste0("IMPORTANCE_consensus_", importance_regime, ".csv")),
            row.names = FALSE)
  cat("\nTop 15 haplotypes by median rank across the", length(unique(long$model)), "retained models:\n")
  print(as.data.frame(head(consensus, 15)), row.names = FALSE, digits = 3)
}

# =============================================================================
# VALIDATION ACCESSIONS MISCLASSIFIED BY EVERY MODEL
# =============================================================================
if ("holdout150" %in% regimes_to_run) {
  wrong_all <- Reduce(`&`, lapply(seeds_to_run, function(s) {
    fr <- read.csv(file.path(out_dir, "holdout150", paste0("seed_", s), "fold_results.csv"),
                   stringsAsFactors = FALSE)
    rowSums(fr[, paste0(model_names, "_correct")]) == 0
  }))
  fr <- read.csv(file.path(out_dir, "holdout150", "seed_1", "fold_results.csv"),
                 stringsAsFactors = FALSE)
  miss <- fr[wrong_all, c("accession", "truth", "has_unseen_haplotype")]
  hap  <- as.data.frame(lapply(df0[miss$accession, haplo_cols], as.character))
  names(hap) <- sub("^Chr ", "", haplo_cols)
  miss <- cbind(miss, hap)
  write.csv(miss, file.path(out_dir, "MISCLASSIFIED_by_all_models_holdout150.csv"), row.names = FALSE)
  cat("\nValidation accessions misclassified by every model in every seed:", nrow(miss), "\n")
  print(miss[, c("accession", "truth", "has_unseen_haplotype")], row.names = FALSE)
}

cat("\nOutputs written under:", normalizePath(out_dir), "\n")
