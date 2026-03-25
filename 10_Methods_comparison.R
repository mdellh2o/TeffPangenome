###############################################################################
# Teff haplotype -> phenotype (White vs Brown)
# FULL “FAIR” NESTED CV + 10 SEEDS + GLOBAL FEATURE-DRIVER SUMMARIES
#
# - OUTER CV: LOOCV
# - INNER CV: Stratified K-fold on outer-train only (tuning)
# - Tuning metric: Balanced accuracy
#
# Encoding:
# - FULL one-hot for ALL factor levels (incl. "0") using identity contrasts
#
# Standardization:
# - For SVM/NB/CART/GBM: fold-wise (train stats only)
# - For glmnet: standardize=FALSE vs TRUE inside glmnet
#
# Models (each with nostd + std):
# - LASSO logistic (glmnet alpha=1)
# - RIDGE logistic (glmnet alpha=0)
# - SVM linear (e1071)
# - SVM polynomial (e1071)
# - Naive Bayes (e1071::naiveBayes)
# - CART (rpart)
# - GBM (gbm)
#
# Extra:
# - Repeat whole experiment for seeds_to_run (default 10)
# - For each seed, compute:
#   (i) comparison table (acc/balacc/sens/spec)
#   (ii) per-model #features "significantly different from 0" across outer folds
#        (heuristic t-test across folds + BH-FDR)
# - Aggregate over seeds:
#   mean/sd of metrics + mean/sd/min/max of n_significant_features
#
# NEW:
# - Extract global feature-level drivers across OUTER FOLDS and SEEDS
# - For signed models: mean sign, nonzero frequency, sign consistency, top-k freq
# - For unsigned models: mean importance, top-k freq, stability
# - Create consensus summary across models
###############################################################################

# ----------------------------
# CONFIG
# ----------------------------
data_file <- "haplo_onQuncho.txt"

# run multiple seeds
seeds_to_run <- 1:10

# Inner CV
innerK_base <- 5

# Hyperparameter grids
svm_cost_grid        <- c(0.1, 1, 10, 100)
svm_poly_cost_grid   <- c(0.1, 1, 10)
svm_poly_degree_grid <- c(2, 3)
svm_poly_coef0_grid  <- c(0, 1)

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

# Importance settings
svm_top_k       <- 10
svm_poly_nperm  <- 10
nb_top_k        <- 10
cart_top_k      <- 10
gbm_top_k       <- 10
global_top_k    <- 10

# significance threshold for "features != 0" heuristic
alpha_sig <- 0.05

# output base folder
out_dir_base <- "loocv_outputs_runs"
dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# PACKAGES
# ----------------------------
pkgs <- c("data.table","dplyr","tidyr","purrr","glmnet","e1071","rpart","gbm")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# =============================================================================
# METRICS + HELPERS
# =============================================================================
clf_metrics <- function(y_true, y_pred, positive_class = "Brown") {
  y_true <- droplevels(factor(y_true))
  y_pred <- factor(y_pred, levels = levels(y_true))
  
  if (nlevels(y_true) < 2) {
    acc <- mean(y_true == y_pred)
    conf <- table(truth = y_true, pred = y_pred)
    return(list(
      acc = acc, bal_acc = NA_real_,
      sensitivity = NA_real_, specificity = NA_real_,
      conf = conf
    ))
  }
  
  neg_class <- setdiff(levels(y_true), positive_class)[1]
  TP <- sum(y_true == positive_class & y_pred == positive_class)
  TN <- sum(y_true == neg_class     & y_pred == neg_class)
  FP <- sum(y_true == neg_class     & y_pred == positive_class)
  FN <- sum(y_true == positive_class & y_pred == neg_class)
  
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  
  list(
    acc = acc,
    bal_acc = bal_acc,
    sensitivity = sens,
    specificity = spec,
    conf = table(truth = y_true, pred = y_pred)
  )
}

bal_acc_safe <- function(y_true, y_pred) {
  m <- clf_metrics(y_true, y_pred, positive_class = "Brown")$bal_acc
  if (is.na(m)) 0 else m
}

standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sigma <- apply(Xtr, 2, sd)
  sigma[is.na(sigma) | sigma == 0] <- 1
  Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/")
  Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sigma, "/")
  list(Xtr = Xtr_s, Xte = Xte_s, mu = mu, sigma = sigma)
}

make_stratified_folds <- function(y, K = 5, seed = 1) {
  set.seed(seed)
  y <- droplevels(factor(y))
  n <- length(y)
  fold_id <- rep(NA_integer_, n)
  
  for (cl in levels(y)) {
    idx <- which(y == cl)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(1:K, length.out = length(idx))
  }
  fold_id
}

prep_std <- function(Xtr, Xva, standardize_flag) {
  if (!standardize_flag) {
    return(list(
      Xtr = Xtr,
      Xva = Xva,
      sigma = rep(1, ncol(Xtr))
    ))
  }
  st <- standardize_train_test(Xtr, Xva)
  list(Xtr = st$Xtr, Xva = st$Xte, sigma = st$sigma)
}

# =============================================================================
# INNER-CV TUNERS (balanced accuracy)
# =============================================================================
tune_glmnet_innerK <- function(Xtr, ytr_bin, ytr_fac, K, standardize_flag,
                               alpha_glmnet, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  
  fit0 <- glmnet::glmnet(
    x = Xtr, y = ytr_bin,
    family = "binomial", alpha = alpha_glmnet,
    standardize = standardize_flag,
    intercept = FALSE
  )
  lambda_grid <- fit0$lambda
  L <- length(lambda_grid)
  
  fold_scores <- matrix(NA_real_, nrow = K, ncol = L)
  
  for (k in 1:K) {
    va <- which(folds == k)
    tr <- which(folds != k)
    if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_bin[tr])) < 2) next
    
    Xk_tr <- Xtr[tr, , drop = FALSE]
    Xk_va <- Xtr[va, , drop = FALSE]
    yk_tr <- ytr_bin[tr]
    yk_va <- ytr_fac[va]
    
    fit <- glmnet::glmnet(
      x = Xk_tr, y = yk_tr,
      family = "binomial", alpha = alpha_glmnet,
      lambda = lambda_grid,
      standardize = standardize_flag,
      intercept = FALSE
    )
    
    pmat <- predict(fit, newx = Xk_va, type = "response")
    used <- fit$lambda
    idx_used <- match(used, lambda_grid)
    
    for (j in seq_along(used)) {
      probs <- as.numeric(pmat[, j])
      pred  <- ifelse(probs >= 0.5, "Brown", "White")
      fold_scores[k, idx_used[j]] <- bal_acc_safe(yk_va, pred)
    }
  }
  
  mean_sc <- colMeans(fold_scores, na.rm = TRUE)
  sd_sc   <- apply(fold_scores, 2, sd, na.rm = TRUE)
  n_eff   <- colSums(!is.na(fold_scores))
  se_sc   <- sd_sc / sqrt(pmax(n_eff, 1))
  
  best_mean <- max(mean_sc, na.rm = TRUE)
  j_best    <- which.max(mean_sc)
  
  thr <- best_mean - se_sc[j_best]
  cand <- which(mean_sc >= thr)
  if (length(cand) == 0) cand <- j_best
  j_1se <- max(cand)
  
  list(
    lambda = lambda_grid[j_1se],
    lambda_grid = lambda_grid,
    mean_scores = mean_sc
  )
}

tune_svm_linear_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                                   cost_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  scores <- setNames(rep(NA_real_, length(cost_grid)), as.character(cost_grid))
  
  for (cst in cost_grid) {
    fold_sc <- rep(NA_real_, K)
    
    for (k in 1:K) {
      va <- which(folds == k)
      tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      
      Xk_tr <- Xtr[tr, , drop = FALSE]
      Xk_va <- Xtr[va, , drop = FALSE]
      yk_tr <- droplevels(ytr_fac[tr])
      yk_va <- ytr_fac[va]
      
      st <- prep_std(Xk_tr, Xk_va, standardize_flag)
      fit <- e1071::svm(
        x = st$Xtr, y = yk_tr,
        kernel = "linear",
        cost = cst,
        scale = FALSE
      )
      pr  <- predict(fit, st$Xva)
      fold_sc[k] <- bal_acc_safe(yk_va, pr)
    }
    scores[as.character(cst)] <- mean(fold_sc, na.rm = TRUE)
  }
  
  best_cost <- as.numeric(names(which.max(scores)))
  list(cost = best_cost, mean_scores = scores)
}

tune_svm_poly_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                                 cost_grid, degree_grid, coef0 = 0, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(cost = cost_grid, degree = degree_grid, stringsAsFactors = FALSE)
  key  <- paste0("c=", grid$cost, "_d=", grid$degree)
  scores <- setNames(rep(NA_real_, nrow(grid)), key)
  
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    
    for (k in 1:K) {
      va <- which(folds == k)
      tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      
      Xk_tr <- Xtr[tr, , drop = FALSE]
      Xk_va <- Xtr[va, , drop = FALSE]
      yk_tr <- droplevels(ytr_fac[tr])
      yk_va <- ytr_fac[va]
      
      st <- prep_std(Xk_tr, Xk_va, standardize_flag)
      fit <- e1071::svm(
        x = st$Xtr, y = yk_tr,
        kernel = "polynomial",
        cost = grid$cost[g],
        degree = grid$degree[g],
        coef0 = coef0,
        scale = FALSE
      )
      pr  <- predict(fit, st$Xva)
      fold_sc[k] <- bal_acc_safe(yk_va, pr)
    }
    scores[key[g]] <- mean(fold_sc, na.rm = TRUE)
  }
  
  best_row <- grid[which.max(scores), , drop = FALSE]
  list(
    cost = best_row$cost,
    degree = best_row$degree,
    coef0 = coef0,
    mean_scores = scores
  )
}

tune_nb_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                           prior_grid, laplace_grid, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(prior_type = prior_grid, laplace = laplace_grid, stringsAsFactors = FALSE)
  key  <- paste0("p=", grid$prior_type, "_l=", grid$laplace)
  scores <- setNames(rep(NA_real_, nrow(grid)), key)
  
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    
    for (k in 1:K) {
      va <- which(folds == k)
      tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      
      Xk_tr <- Xtr[tr, , drop = FALSE]
      Xk_va <- Xtr[va, , drop = FALSE]
      yk_tr <- droplevels(ytr_fac[tr])
      yk_va <- ytr_fac[va]
      
      st <- prep_std(Xk_tr, Xk_va, standardize_flag)
      
      if (grid$prior_type[g] == "uniform") {
        pr0 <- rep(1 / nlevels(yk_tr), nlevels(yk_tr))
        names(pr0) <- levels(yk_tr)
        fit <- e1071::naiveBayes(
          x = st$Xtr, y = yk_tr,
          prior = pr0,
          laplace = grid$laplace[g]
        )
      } else {
        fit <- e1071::naiveBayes(
          x = st$Xtr, y = yk_tr,
          laplace = grid$laplace[g]
        )
      }
      
      pr <- predict(fit, st$Xva)
      fold_sc[k] <- bal_acc_safe(yk_va, pr)
    }
    scores[key[g]] <- mean(fold_sc, na.rm = TRUE)
  }
  
  best_idx <- which.max(scores)
  list(
    prior_type = grid$prior_type[best_idx],
    laplace = grid$laplace[best_idx],
    mean_scores = scores
  )
}

tune_cart_innerK <- function(Xtr, ytr_fac, K, standardize_flag,
                             cp_grid, maxdepth, minsplit, seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  scores <- setNames(rep(NA_real_, length(cp_grid)), as.character(cp_grid))
  
  for (cpv in cp_grid) {
    fold_sc <- rep(NA_real_, K)
    
    for (k in 1:K) {
      va <- which(folds == k)
      tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      
      Xk_tr <- Xtr[tr, , drop = FALSE]
      Xk_va <- Xtr[va, , drop = FALSE]
      yk_tr <- droplevels(ytr_fac[tr])
      yk_va <- ytr_fac[va]
      
      st <- prep_std(Xk_tr, Xk_va, standardize_flag)
      
      df_tr <- data.frame(y = yk_tr, as.data.frame(st$Xtr))
      ctrl <- rpart::rpart.control(
        cp = cpv,
        maxdepth = maxdepth,
        minsplit = minsplit,
        xval = 0
      )
      fit <- rpart::rpart(
        y ~ ., data = df_tr,
        method = "class",
        control = ctrl
      )
      
      pr <- predict(fit, newdata = as.data.frame(st$Xva), type = "class")
      fold_sc[k] <- bal_acc_safe(yk_va, pr)
    }
    scores[as.character(cpv)] <- mean(fold_sc, na.rm = TRUE)
  }
  
  best_cp <- as.numeric(names(which.max(scores)))
  list(cp = best_cp, mean_scores = scores)
}

tune_gbm_innerK <- function(Xtr, ytr_bin, ytr_fac, K, standardize_flag,
                            depth_grid, trees_grid, shrinkage, minobs, bagfrac,
                            seed = 1) {
  folds <- make_stratified_folds(ytr_fac, K = K, seed = seed)
  grid <- expand.grid(depth = depth_grid, trees = trees_grid, stringsAsFactors = FALSE)
  key  <- paste0("d=", grid$depth, "_t=", grid$trees)
  scores <- setNames(rep(NA_real_, nrow(grid)), key)
  
  max_t <- max(trees_grid)
  
  for (g in seq_len(nrow(grid))) {
    fold_sc <- rep(NA_real_, K)
    dep <- grid$depth[g]
    tt  <- grid$trees[g]
    
    for (k in 1:K) {
      va <- which(folds == k)
      tr <- which(folds != k)
      if (length(unique(ytr_fac[va])) < 2 || length(unique(ytr_fac[tr])) < 2) next
      
      Xk_tr <- Xtr[tr, , drop = FALSE]
      Xk_va <- Xtr[va, , drop = FALSE]
      yk_trb <- ytr_bin[tr]
      yk_vaf <- ytr_fac[va]
      
      st <- prep_std(Xk_tr, Xk_va, standardize_flag)
      df_tr <- data.frame(y = yk_trb, as.data.frame(st$Xtr))
      
      set.seed(seed + 10000 * g + 100 * k + if (standardize_flag) 1 else 0)
      
      fit <- gbm::gbm(
        y ~ .,
        data = df_tr,
        distribution = "bernoulli",
        n.trees = max_t,
        interaction.depth = dep,
        shrinkage = shrinkage,
        n.minobsinnode = minobs,
        bag.fraction = bagfrac,
        train.fraction = 1.0,
        cv.folds = 0,
        verbose = FALSE
      )
      
      p_va <- predict(fit, newdata = as.data.frame(st$Xva), n.trees = tt, type = "response")
      pr   <- ifelse(p_va >= 0.5, "Brown", "White")
      fold_sc[k] <- bal_acc_safe(yk_vaf, pr)
    }
    
    scores[key[g]] <- mean(fold_sc, na.rm = TRUE)
  }
  
  best_idx <- which.max(scores)
  list(
    depth = grid$depth[best_idx],
    trees = grid$trees[best_idx],
    mean_scores = scores
  )
}

# =============================================================================
# IMPORTANCE EXTRACTORS
# =============================================================================
importance_glmnet <- function(fit, lambda, full_feat_names) {
  b <- as.matrix(coef(fit, s = lambda))
  rn <- rownames(b)
  if ("(Intercept)" %in% rn) b <- b[setdiff(rn, "(Intercept)"), , drop = FALSE]
  v <- rep(0, length(full_feat_names))
  names(v) <- full_feat_names
  common <- intersect(rownames(b), full_feat_names)
  v[common] <- b[common, 1]
  v
}

importance_svm_linear <- function(svm_fit, Xtr_fit, ytr,
                                  positive_class = "Brown", sigma) {
  w_scaled <- drop(t(svm_fit$coefs) %*% svm_fit$SV)
  names(w_scaled) <- colnames(Xtr_fit)
  
  pr_tr <- predict(svm_fit, Xtr_fit, decision.values = TRUE)
  dv_tr <- as.numeric(attr(pr_tr, "decision.values"))
  mean_pos <- mean(dv_tr[ytr == positive_class], na.rm = TRUE)
  mean_neg <- mean(dv_tr[ytr != positive_class], na.rm = TRUE)
  sgn <- ifelse((mean_pos - mean_neg) >= 0, 1, -1)
  
  (w_scaled / sigma) * sgn
}

importance_svm_poly_perm <- function(fit, Xtr_fit, ytr, nperm = 10, seed = 1) {
  set.seed(seed)
  base_pr  <- predict(fit, Xtr_fit)
  base_bal <- bal_acc_safe(ytr, base_pr)
  
  p <- ncol(Xtr_fit)
  imp <- rep(0, p)
  names(imp) <- colnames(Xtr_fit)
  
  for (j in seq_len(p)) {
    drops <- numeric(nperm)
    for (b in 1:nperm) {
      Xp <- Xtr_fit
      Xp[, j] <- sample(Xp[, j], replace = FALSE)
      prp <- predict(fit, Xp)
      drops[b] <- base_bal - bal_acc_safe(ytr, prp)
    }
    imp[j] <- mean(pmax(drops, 0))
  }
  imp
}

importance_nb_effect <- function(Xtr_fit, ytr, positive_class = "Brown") {
  Xpos <- Xtr_fit[ytr == positive_class, , drop = FALSE]
  Xneg <- Xtr_fit[ytr != positive_class, , drop = FALSE]
  mu_pos <- colMeans(Xpos)
  mu_neg <- colMeans(Xneg)
  sd_pos <- apply(Xpos, 2, sd)
  sd_neg <- apply(Xneg, 2, sd)
  sd_pos[is.na(sd_pos)] <- 0
  sd_neg[is.na(sd_neg)] <- 0
  sd_pool <- sqrt((sd_pos^2 + sd_neg^2) / 2)
  abs(mu_pos - mu_neg) / (sd_pool + 1e-8)
}

importance_cart_vi <- function(fit, feat_names) {
  vi <- fit$variable.importance
  v <- rep(0, length(feat_names))
  names(v) <- feat_names
  if (!is.null(vi) && length(vi) > 0) {
    common <- intersect(names(vi), feat_names)
    v[common] <- as.numeric(vi[common])
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

importance_gbm_rel <- function(fit, n.trees, feat_names) {
  rel <- suppressWarnings(summary(fit, plotit = FALSE, n.trees = n.trees))
  v <- rep(0, length(feat_names))
  names(v) <- feat_names
  if (!is.null(rel) && nrow(rel) > 0) {
    common <- intersect(as.character(rel$var), feat_names)
    v[common] <- rel$rel.inf[match(common, as.character(rel$var))]
    if (sum(v) > 0) v <- v / sum(v)
  }
  v
}

# =============================================================================
# IMPORTANCE SUMMARIES
# =============================================================================
summarise_signed_importance <- function(M, top_k = 10) {
  mean_w <- colMeans(M, na.rm = TRUE)
  sd_w   <- apply(M, 2, sd, na.rm = TRUE)
  stab   <- abs(mean_w) / (sd_w + 1e-8)
  
  topk_mat <- matrix(FALSE, nrow(M), ncol(M))
  colnames(topk_mat) <- colnames(M)
  for (i in 1:nrow(M)) {
    wi <- M[i, ]
    ord <- order(abs(wi), decreasing = TRUE)
    k <- min(top_k, length(ord))
    topk_mat[i, ord[1:k]] <- TRUE
  }
  topk_freq <- colMeans(topk_mat)
  
  mean_sign <- sign(mean_w)
  sign_cons <- sapply(seq_along(mean_w), function(j) {
    if (mean_sign[j] == 0) return(NA_real_)
    mean(sign(M[, j]) == mean_sign[j], na.rm = TRUE)
  })
  
  data.frame(
    term = colnames(M),
    mean = as.numeric(mean_w),
    sd = as.numeric(sd_w),
    mean_abs = abs(as.numeric(mean_w)),
    sign_consistency = as.numeric(sign_cons),
    topk_freq = as.numeric(topk_freq),
    stability_score = as.numeric(stab),
    direction = ifelse(mean_w > 0, "Brown", "White"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(desc(stability_score), desc(mean_abs))
}

summarise_unsigned_importance <- function(M, top_k = 10) {
  mean_i <- colMeans(M, na.rm = TRUE)
  sd_i   <- apply(M, 2, sd, na.rm = TRUE)
  stab   <- mean_i / (sd_i + 1e-8)
  
  topk_mat <- matrix(FALSE, nrow(M), ncol(M))
  colnames(topk_mat) <- colnames(M)
  for (i in 1:nrow(M)) {
    vi <- M[i, ]
    ord <- order(vi, decreasing = TRUE)
    k <- min(top_k, length(ord))
    topk_mat[i, ord[1:k]] <- TRUE
  }
  topk_freq <- colMeans(topk_mat)
  
  data.frame(
    term = colnames(M),
    mean = as.numeric(mean_i),
    sd = as.numeric(sd_i),
    stability_score = as.numeric(stab),
    topk_freq = as.numeric(topk_freq),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(desc(stability_score), desc(mean))
}

# =============================================================================
# “SIGNIFICANT FEATURES” COUNT (heuristic across outer folds)
# =============================================================================
count_sig_features <- function(M, signed = TRUE, alpha = 0.05, eps = 1e-12) {
  xbar <- colMeans(M, na.rm = TRUE)
  sdev <- apply(M, 2, sd, na.rm = TRUE)
  n_eff <- colSums(!is.na(M))
  se <- sdev / sqrt(pmax(n_eff, 1))
  tstat <- xbar / (se + eps)
  df <- pmax(n_eff - 1, 1)
  
  if (signed) {
    p <- 2 * pt(-abs(tstat), df = df)
  } else {
    p <- 1 - pt(tstat, df = df)
  }
  q <- p.adjust(p, method = "BH")
  sum(q < alpha, na.rm = TRUE)
}

# =============================================================================
# FEATURE USAGE SUMMARY
# =============================================================================
add_raw_names <- function(df, safe_to_raw, term_col = "term") {
  df$raw_term <- safe_to_raw[df[[term_col]]]
  df$raw_term[is.na(df$raw_term)] <- df[[term_col]][is.na(df$raw_term)]
  df
}

summarise_sparse_usage <- function(M, tau = 0.8, eps = 1e-12) {
  nz <- abs(M) > eps
  freq <- colMeans(nz, na.rm = TRUE)
  nnz_per_fold <- rowSums(nz, na.rm = TRUE)
  
  sign_cons <- rep(NA_real_, ncol(M))
  for (j in seq_len(ncol(M))) {
    idx <- which(nz[, j])
    if (length(idx) >= 2) {
      s <- sign(M[idx, j])
      s0 <- sign(mean(M[idx, j]))
      sign_cons[j] <- mean(s == s0)
    }
  }
  
  list(
    per_fold = data.frame(
      fold = rownames(M),
      n_nonzero = nnz_per_fold,
      stringsAsFactors = FALSE
    ),
    by_feature = data.frame(
      term = colnames(M),
      sel_freq = freq,
      sign_consistency = sign_cons,
      stringsAsFactors = FALSE
    ),
    n_consistent = sum(freq >= tau, na.rm = TRUE)
  )
}

effective_feature_count <- function(M, eps = 1e-12) {
  apply(M, 1, function(w) {
    w <- as.numeric(w)
    l1 <- sum(abs(w), na.rm = TRUE)
    l2 <- sqrt(sum(w^2, na.rm = TRUE))
    if (l2 < eps) return(0)
    (l1^2) / (l2^2 + eps)
  })
}

approx_significance_across_folds <- function(M, signed = TRUE, alpha = 0.05,
                                             eps = 1e-12) {
  xbar <- colMeans(M, na.rm = TRUE)
  sdev <- apply(M, 2, sd, na.rm = TRUE)
  n_eff <- colSums(!is.na(M))
  se <- sdev / sqrt(pmax(n_eff, 1))
  tstat <- xbar / (se + eps)
  df <- pmax(n_eff - 1, 1)
  
  if (signed) {
    p <- 2 * pt(-abs(tstat), df = df)
  } else {
    p <- 1 - pt(tstat, df = df)
  }
  q <- p.adjust(p, method = "BH")
  
  out <- data.frame(
    term = colnames(M),
    mean = as.numeric(xbar),
    sd = as.numeric(sdev),
    n_eff = as.integer(n_eff),
    t = as.numeric(tstat),
    p = as.numeric(p),
    q = as.numeric(q),
    significant = q < alpha,
    stringsAsFactors = FALSE
  )
  out[order(out$q, -abs(out$mean)), , drop = FALSE]
}

build_usage_summary <- function(IMP, name_map, out_dir, tau = 0.8, alpha = 0.05) {
  safe_to_raw <- setNames(name_map$raw, name_map$safe)
  
  sparse_methods <- c("lasso_nostd", "lasso_std")
  signed_dense_methods <- c("svm_lin_nostd", "svm_lin_std", "l2log_nostd", "l2log_std")
  unsigned_methods <- c("svm_poly_nostd","svm_poly_std",
                        "nb_nostd","nb_std",
                        "cart_nostd","cart_std",
                        "gbm_nostd","gbm_std")
  
  usage_summary <- list()
  
  for (m in sparse_methods) {
    M <- IMP[[m]]
    su <- summarise_sparse_usage(M, tau = tau)
    
    sig <- approx_significance_across_folds(M, signed = TRUE, alpha = alpha)
    sig <- add_raw_names(sig, safe_to_raw, "term")
    
    usage_summary[[m]] <- data.frame(
      model = m,
      mean_nonzero_per_fold = mean(su$per_fold$n_nonzero),
      median_nonzero_per_fold = median(su$per_fold$n_nonzero),
      consistent_features_freq_ge_tau = su$n_consistent,
      approx_significant_q_lt_alpha = sum(sig$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    write.csv(add_raw_names(su$by_feature, safe_to_raw, "term"),
              file.path(out_dir, paste0(m, "_selection_frequency.csv")),
              row.names = FALSE)
    write.csv(sig,
              file.path(out_dir, paste0(m, "_approx_significance.csv")),
              row.names = FALSE)
  }
  
  for (m in signed_dense_methods) {
    M <- IMP[[m]]
    neff <- effective_feature_count(M)
    
    sig <- approx_significance_across_folds(M, signed = TRUE, alpha = alpha)
    sig <- add_raw_names(sig, safe_to_raw, "term")
    
    usage_summary[[m]] <- data.frame(
      model = m,
      mean_effective_features_Neff = mean(neff),
      median_effective_features_Neff = median(neff),
      approx_significant_q_lt_alpha = sum(sig$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    write.csv(data.frame(fold = rownames(M), Neff = neff),
              file.path(out_dir, paste0(m, "_effective_feature_count.csv")),
              row.names = FALSE)
    write.csv(sig,
              file.path(out_dir, paste0(m, "_approx_significance.csv")),
              row.names = FALSE)
  }
  
  for (m in unsigned_methods) {
    M <- IMP[[m]]
    neff <- effective_feature_count(M)
    
    sig <- approx_significance_across_folds(M, signed = FALSE, alpha = alpha)
    sig <- add_raw_names(sig, safe_to_raw, "term")
    
    usage_summary[[m]] <- data.frame(
      model = m,
      mean_effective_features_Neff = mean(neff),
      median_effective_features_Neff = median(neff),
      approx_significant_q_lt_alpha = sum(sig$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    write.csv(data.frame(fold = rownames(M), Neff = neff),
              file.path(out_dir, paste0(m, "_effective_feature_count.csv")),
              row.names = FALSE)
    write.csv(sig,
              file.path(out_dir, paste0(m, "_approx_significance.csv")),
              row.names = FALSE)
  }
  
  usage_summary_df <- data.table::rbindlist(usage_summary, fill = TRUE)
  usage_summary_df <- as.data.frame(usage_summary_df)
  usage_summary_df <- usage_summary_df[order(usage_summary_df$model), ]
  write.csv(usage_summary_df, file.path(out_dir, "feature_usage_summary.csv"), row.names = FALSE)
  usage_summary_df
}

# =============================================================================
# GLOBAL FEATURE-DRIVER HELPERS
# =============================================================================
pretty_feature_name <- function(x) {
  x2 <- x
  x2 <- gsub("^Chr\\s*([0-9]+[A-Z])([0-9]+)$", "Chr \\1 = \\2", x2)
  x2 <- gsub("^Chr\\.([0-9]+[A-Z])\\.([0-9]+)$", "Chr \\1 = \\2", x2)
  x2 <- gsub("^Chr\\s+([0-9]+[A-Z])\\s+([0-9]+)$", "Chr \\1 = \\2", x2)
  x2
}

collect_importance_across_seeds <- function(out_dir_base, seeds_to_run, model_name) {
  mats <- list()
  for (s in seeds_to_run) {
    f <- file.path(out_dir_base, paste0("seed_", s), "loocv_all_objects.rds")
    obj <- readRDS(f)
    M <- obj$importance_matrices[[model_name]]
    rownames(M) <- paste0("seed", s, "__", rownames(M))
    mats[[as.character(s)]] <- M
  }
  do.call(rbind, mats)
}

build_global_feature_driver_summaries <- function(out_dir_base, seeds_to_run,
                                                  top_k = 10) {
  all_models <- c(
    "lasso_nostd","lasso_std",
    "l2log_nostd","l2log_std",
    "svm_lin_nostd","svm_lin_std",
    "svm_poly_nostd","svm_poly_std",
    "nb_nostd","nb_std",
    "cart_nostd","cart_std",
    "gbm_nostd","gbm_std"
  )
  
  signed_models <- c(
    "lasso_nostd","lasso_std",
    "l2log_nostd","l2log_std",
    "svm_lin_nostd","svm_lin_std"
  )
  unsigned_models <- setdiff(all_models, signed_models)
  
  name_map <- read.csv(
    file.path(out_dir_base, "seed_1", "feature_name_map.csv"),
    stringsAsFactors = FALSE
  )
  safe_to_raw <- setNames(name_map$raw, name_map$safe)
  
  global_feature_summaries <- list()
  
  for (m in signed_models) {
    M <- collect_importance_across_seeds(out_dir_base, seeds_to_run, m)
    
    mean_w   <- colMeans(M, na.rm = TRUE)
    sd_w     <- apply(M, 2, sd, na.rm = TRUE)
    mean_abs <- colMeans(abs(M), na.rm = TRUE)
    nz_freq  <- colMeans(abs(M) > 1e-12, na.rm = TRUE)
    
    topk_mat <- matrix(FALSE, nrow = nrow(M), ncol = ncol(M))
    colnames(topk_mat) <- colnames(M)
    for (i in 1:nrow(M)) {
      ord <- order(abs(M[i, ]), decreasing = TRUE)
      k <- min(top_k, length(ord))
      topk_mat[i, ord[1:k]] <- TRUE
    }
    topk_freq <- colMeans(topk_mat, na.rm = TRUE)
    
    sign_cons <- sapply(seq_along(mean_w), function(j) {
      vals <- M[, j]
      vals <- vals[abs(vals) > 1e-12]
      if (length(vals) < 2) return(NA_real_)
      s0 <- sign(mean(vals))
      mean(sign(vals) == s0, na.rm = TRUE)
    })
    
    stab <- abs(mean_w) / (sd_w + 1e-8)
    
    df_sum <- data.frame(
      model = m,
      term_safe = colnames(M),
      term_raw = safe_to_raw[colnames(M)],
      mean = as.numeric(mean_w),
      mean_abs = as.numeric(mean_abs),
      sd = as.numeric(sd_w),
      nonzero_freq = as.numeric(nz_freq),
      topk_freq = as.numeric(topk_freq),
      sign_consistency = as.numeric(sign_cons),
      stability_score = as.numeric(stab),
      direction = ifelse(mean_w > 0, "Brown", ifelse(mean_w < 0, "White", "None")),
      stringsAsFactors = FALSE
    )
    
    df_sum$term_raw[is.na(df_sum$term_raw)] <- df_sum$term_safe[is.na(df_sum$term_raw)]
    df_sum$feature_label <- pretty_feature_name(df_sum$term_raw)
    
    df_sum <- df_sum[order(-df_sum$stability_score, -df_sum$mean_abs), ]
    global_feature_summaries[[m]] <- df_sum
    
    write.csv(
      df_sum,
      file.path(out_dir_base, paste0("GLOBAL_feature_drivers_", m, ".csv")),
      row.names = FALSE
    )
  }
  
  for (m in unsigned_models) {
    M <- collect_importance_across_seeds(out_dir_base, seeds_to_run, m)
    
    mean_i  <- colMeans(M, na.rm = TRUE)
    sd_i    <- apply(M, 2, sd, na.rm = TRUE)
    nz_freq <- colMeans(M > 1e-12, na.rm = TRUE)
    
    topk_mat <- matrix(FALSE, nrow = nrow(M), ncol = ncol(M))
    colnames(topk_mat) <- colnames(M)
    for (i in 1:nrow(M)) {
      ord <- order(M[i, ], decreasing = TRUE)
      k <- min(top_k, length(ord))
      topk_mat[i, ord[1:k]] <- TRUE
    }
    topk_freq <- colMeans(topk_mat, na.rm = TRUE)
    
    stab <- mean_i / (sd_i + 1e-8)
    
    df_sum <- data.frame(
      model = m,
      term_safe = colnames(M),
      term_raw = safe_to_raw[colnames(M)],
      mean_importance = as.numeric(mean_i),
      sd = as.numeric(sd_i),
      nonzero_freq = as.numeric(nz_freq),
      topk_freq = as.numeric(topk_freq),
      stability_score = as.numeric(stab),
      stringsAsFactors = FALSE
    )
    
    df_sum$term_raw[is.na(df_sum$term_raw)] <- df_sum$term_safe[is.na(df_sum$term_raw)]
    df_sum$feature_label <- pretty_feature_name(df_sum$term_raw)
    
    df_sum <- df_sum[order(-df_sum$stability_score, -df_sum$mean_importance), ]
    global_feature_summaries[[m]] <- df_sum
    
    write.csv(
      df_sum,
      file.path(out_dir_base, paste0("GLOBAL_feature_drivers_", m, ".csv")),
      row.names = FALSE
    )
  }
  
  consensus_list <- list()
  
  for (m in names(global_feature_summaries)) {
    dfm <- global_feature_summaries[[m]]
    
    if (m %in% signed_models) {
      tmp <- dfm[, c("feature_label","topk_freq","nonzero_freq","stability_score","direction")]
      names(tmp) <- c("feature_label","topk_freq","nonzero_freq","stability_score","direction")
    } else {
      tmp <- dfm[, c("feature_label","topk_freq","nonzero_freq","stability_score")]
      tmp$direction <- NA_character_
    }
    
    tmp$model <- m
    consensus_list[[m]] <- tmp
  }
  
  consensus_long <- dplyr::bind_rows(consensus_list)
  
  consensus_summary <- consensus_long %>%
    dplyr::group_by(feature_label) %>%
    dplyr::summarise(
      mean_topk_freq = mean(topk_freq, na.rm = TRUE),
      mean_nonzero_freq = mean(nonzero_freq, na.rm = TRUE),
      mean_stability = mean(stability_score, na.rm = TRUE),
      n_models_appearing = sum(topk_freq > 0, na.rm = TRUE),
      brown_votes = sum(direction == "Brown", na.rm = TRUE),
      white_votes = sum(direction == "White", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      consensus_direction = dplyr::case_when(
        brown_votes > white_votes ~ "Brown",
        white_votes > brown_votes ~ "White",
        TRUE ~ "Mixed/NA"
      )
    ) %>%
    dplyr::arrange(
      dplyr::desc(mean_topk_freq),
      dplyr::desc(n_models_appearing),
      dplyr::desc(mean_stability)
    )
  
  write.csv(
    consensus_long,
    file.path(out_dir_base, "GLOBAL_feature_consensus_long.csv"),
    row.names = FALSE
  )
  write.csv(
    consensus_summary,
    file.path(out_dir_base, "GLOBAL_feature_consensus_summary.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      global_feature_summaries = global_feature_summaries,
      consensus_long = consensus_long,
      consensus_summary = consensus_summary
    ),
    file.path(out_dir_base, "GLOBAL_feature_driver_objects.rds")
  )
  
  invisible(list(
    global_feature_summaries = global_feature_summaries,
    consensus_long = consensus_long,
    consensus_summary = consensus_summary
  ))
}

# =============================================================================
# ONE COMPLETE RUN (ONE SEED)
# =============================================================================
run_one_seed <- function(seed, out_dir, data_file) {
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------
  # LOAD DATA
  # ----------------------------
  onquncho <- data.table::fread(data_file, data.table = FALSE)
  stopifnot(ncol(onquncho) >= 9)
  
  colnames(onquncho)[3:9] <- c("Chr 1B","Chr 3B","Chr 4A","Chr 4B","Chr 5A","Chr 6A","Chr 9A")
  rownames(onquncho) <- onquncho[, 1]
  onquncho <- onquncho[, -1, drop = FALSE]
  
  df0 <- onquncho %>%
    dplyr::filter(phenotype %in% c("White","Brown")) %>%
    dplyr::mutate(phenotype = factor(phenotype, levels = c("White","Brown"))) %>%
    dplyr::mutate(dplyr::across(-phenotype, ~ factor(as.character(.x))))
  
  haplo_cols <- setdiff(names(df0), "phenotype")
  
  # ----------------------------
  # DESIGN MATRIX — FULL one-hot (includes level "0")
  # ----------------------------
  markers_df <- df0[, haplo_cols, drop = FALSE]
  
  contr_list <- lapply(markers_df, function(z) {
    z <- droplevels(z)
    mm <- diag(nlevels(z))
    rownames(mm) <- levels(z)
    colnames(mm) <- levels(z)
    mm
  })
  
  X_all <- model.matrix(~ . - 1, data = markers_df, contrasts.arg = contr_list)
  
  colnames_raw  <- colnames(X_all)
  colnames_safe <- make.names(colnames_raw, unique = TRUE)
  colnames_safe <- sub("^(.*?)([0-9]+)$", "\\1.\\2", colnames_safe)
  colnames(X_all) <- colnames_safe
  
  name_map <- data.frame(raw = colnames_raw, safe = colnames(X_all), stringsAsFactors = FALSE)
  write.csv(name_map, file.path(out_dir, "feature_name_map.csv"), row.names = FALSE)
  
  y_fac <- df0$phenotype
  y_bin <- ifelse(y_fac == "Brown", 1, 0)
  
  ids <- rownames(df0)
  feat_names <- colnames(X_all)
  
  # ----------------------------
  # OUTER NESTED LOOCV
  # ----------------------------
  n <- nrow(X_all)
  
  pred_store <- list(
    lasso_nostd = rep(NA_character_, n),
    lasso_std   = rep(NA_character_, n),
    l2log_nostd = rep(NA_character_, n),
    l2log_std   = rep(NA_character_, n),
    svm_lin_nostd = rep(NA_character_, n),
    svm_lin_std   = rep(NA_character_, n),
    svm_poly_nostd = rep(NA_character_, n),
    svm_poly_std   = rep(NA_character_, n),
    nb_nostd = rep(NA_character_, n),
    nb_std   = rep(NA_character_, n),
    cart_nostd = rep(NA_character_, n),
    cart_std   = rep(NA_character_, n),
    gbm_nostd = rep(NA_character_, n),
    gbm_std   = rep(NA_character_, n)
  )
  
  hp_store <- lapply(pred_store, function(x) vector("list", n))
  
  IMP <- lapply(pred_store, function(x) {
    matrix(0, nrow = n, ncol = length(feat_names), dimnames = list(ids, feat_names))
  })
  
  for (i in 1:n) {
    test_idx  <- i
    train_idx <- setdiff(seq_len(n), i)
    
    Xtr0 <- X_all[train_idx, , drop = FALSE]
    Xte0 <- X_all[test_idx, , drop = FALSE]
    ytr_fac <- y_fac[train_idx]
    ytr_bin <- y_bin[train_idx]
    
    min_class_n <- min(table(ytr_fac))
    innerK <- min(innerK_base, min_class_n)
    innerK <- max(innerK, 2)
    
    sd_outer <- apply(Xtr0, 2, sd)
    sd_outer[is.na(sd_outer)] <- 0
    keep_outer <- sd_outer > 0
    
    if (!any(keep_outer)) {
      maj <- names(which.max(table(ytr_fac)))
      for (nm in names(pred_store)) pred_store[[nm]][i] <- maj
      next
    }
    
    Xtr <- Xtr0[, keep_outer, drop = FALSE]
    Xte <- Xte0[, keep_outer, drop = FALSE]
    feat_outer <- colnames(Xtr)
    
    # --------------------------
    # LASSO
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "lasso_std" else "lasso_nostd"
      
      tune <- tune_glmnet_innerK(
        Xtr, ytr_bin, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        alpha_glmnet = 1,
        seed = seed + 1000 * i + if (std_flag) 1 else 0
      )
      lambda_best <- tune$lambda
      hp_store[[tag]][[i]] <- list(lambda = lambda_best, innerK = innerK)
      
      fit <- glmnet::glmnet(
        x = Xtr, y = ytr_bin,
        family = "binomial",
        alpha = 1,
        lambda = lambda_best,
        standardize = std_flag,
        intercept = FALSE
      )
      
      p_te <- as.numeric(predict(fit, newx = Xte, type = "response"))
      pred_store[[tag]][i] <- ifelse(p_te >= 0.5, "Brown", "White")
      
      v_keep <- importance_glmnet(fit, lambda_best, full_feat_names = feat_outer)
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- v_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # RIDGE / L2 LOGISTIC
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "l2log_std" else "l2log_nostd"
      
      tune <- tune_glmnet_innerK(
        Xtr, ytr_bin, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        alpha_glmnet = 0,
        seed = seed + 1500 * i + if (std_flag) 1 else 0
      )
      lambda_best <- tune$lambda
      hp_store[[tag]][[i]] <- list(lambda = lambda_best, innerK = innerK)
      
      fit <- glmnet::glmnet(
        x = Xtr, y = ytr_bin,
        family = "binomial",
        alpha = 0,
        lambda = lambda_best,
        standardize = std_flag,
        intercept = FALSE
      )
      
      p_te <- as.numeric(predict(fit, newx = Xte, type = "response"))
      pred_store[[tag]][i] <- ifelse(p_te >= 0.5, "Brown", "White")
      
      v_keep <- importance_glmnet(fit, lambda_best, full_feat_names = feat_outer)
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- v_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # SVM linear
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "svm_lin_std" else "svm_lin_nostd"
      
      tune <- tune_svm_linear_innerK(
        Xtr, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        cost_grid = svm_cost_grid,
        seed = seed + 2000 * i + if (std_flag) 1 else 0
      )
      cost_best <- tune$cost
      hp_store[[tag]][[i]] <- list(cost = cost_best, innerK = innerK)
      
      st_outer <- prep_std(Xtr, Xte, std_flag)
      ytr_drop <- droplevels(ytr_fac)
      
      fit <- e1071::svm(
        x = st_outer$Xtr, y = ytr_drop,
        kernel = "linear",
        cost = cost_best,
        scale = FALSE,
        decision.values = TRUE
      )
      
      pr_te <- predict(fit, st_outer$Xva, decision.values = TRUE)
      pred_store[[tag]][i] <- as.character(pr_te)
      
      sigma <- if (std_flag) st_outer$sigma else rep(1, ncol(st_outer$Xtr))
      w_keep <- importance_svm_linear(
        fit, st_outer$Xtr, ytr_drop,
        "Brown",
        sigma = sigma
      )
      
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- w_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # SVM polynomial
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "svm_poly_std" else "svm_poly_nostd"
      
      tune <- tune_svm_poly_innerK(
        Xtr, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        cost_grid = svm_poly_cost_grid,
        degree_grid = svm_poly_degree_grid,
        coef0 = svm_poly_coef0_grid,
        seed = seed + 3000 * i + if (std_flag) 1 else 0
      )
      hp_store[[tag]][[i]] <- list(
        cost = tune$cost,
        degree = tune$degree,
        coef0 = tune$coef0,
        innerK = innerK
      )
      
      st_outer <- prep_std(Xtr, Xte, std_flag)
      ytr_drop <- droplevels(ytr_fac)
      
      fit <- e1071::svm(
        x = st_outer$Xtr, y = ytr_drop,
        kernel = "polynomial",
        cost = tune$cost,
        degree = tune$degree,
        coef0 = tune$coef0,
        scale = FALSE
      )
      
      pr_te <- predict(fit, st_outer$Xva)
      pred_store[[tag]][i] <- as.character(pr_te)
      
      imp_keep <- importance_svm_poly_perm(
        fit, st_outer$Xtr, ytr_drop,
        nperm = svm_poly_nperm,
        seed = seed + 3100 * i + if (std_flag) 1 else 0
      )
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- imp_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # Naive Bayes
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "nb_std" else "nb_nostd"
      
      tune <- tune_nb_innerK(
        Xtr, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        prior_grid = nb_prior_grid,
        laplace_grid = nb_laplace_grid,
        seed = seed + 4000 * i + if (std_flag) 1 else 0
      )
      hp_store[[tag]][[i]] <- list(
        prior_type = tune$prior_type,
        laplace = tune$laplace,
        innerK = innerK
      )
      
      st_outer <- prep_std(Xtr, Xte, std_flag)
      ytr_drop <- droplevels(ytr_fac)
      
      if (tune$prior_type == "uniform") {
        pr0 <- rep(1 / nlevels(ytr_drop), nlevels(ytr_drop))
        names(pr0) <- levels(ytr_drop)
        fit <- e1071::naiveBayes(
          x = st_outer$Xtr, y = ytr_drop,
          prior = pr0,
          laplace = tune$laplace
        )
      } else {
        fit <- e1071::naiveBayes(
          x = st_outer$Xtr, y = ytr_drop,
          laplace = tune$laplace
        )
      }
      
      pr_te <- predict(fit, st_outer$Xva)
      pred_store[[tag]][i] <- as.character(pr_te)
      
      eff_keep <- importance_nb_effect(st_outer$Xtr, ytr_drop, positive_class = "Brown")
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- eff_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # CART
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "cart_std" else "cart_nostd"
      
      tune <- tune_cart_innerK(
        Xtr, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        cp_grid = cart_cp_grid,
        maxdepth = cart_maxdepth,
        minsplit = cart_minsplit,
        seed = seed + 5000 * i + if (std_flag) 1 else 0
      )
      hp_store[[tag]][[i]] <- list(
        cp = tune$cp,
        maxdepth = cart_maxdepth,
        minsplit = cart_minsplit,
        innerK = innerK
      )
      
      st_outer <- prep_std(Xtr, Xte, std_flag)
      ytr_drop <- droplevels(ytr_fac)
      
      df_tr <- data.frame(y = ytr_drop, as.data.frame(st_outer$Xtr))
      ctrl  <- rpart::rpart.control(
        cp = tune$cp,
        maxdepth = cart_maxdepth,
        minsplit = cart_minsplit,
        xval = 0
      )
      fit <- rpart::rpart(
        y ~ ., data = df_tr,
        method = "class",
        control = ctrl
      )
      
      pr_te <- predict(fit, newdata = as.data.frame(st_outer$Xva), type = "class")
      pred_store[[tag]][i] <- as.character(pr_te)
      
      vi_keep <- importance_cart_vi(fit, feat_outer)
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- vi_keep
      IMP[[tag]][i, ] <- v_full
    }
    
    # --------------------------
    # GBM
    # --------------------------
    for (std_flag in c(FALSE, TRUE)) {
      tag <- if (std_flag) "gbm_std" else "gbm_nostd"
      
      tune <- tune_gbm_innerK(
        Xtr, ytr_bin, ytr_fac,
        K = innerK,
        standardize_flag = std_flag,
        depth_grid = gbm_depth_grid,
        trees_grid = gbm_trees_grid,
        shrinkage = gbm_shrinkage,
        minobs = gbm_minobs,
        bagfrac = gbm_bagfrac,
        seed = seed + 6000 * i + if (std_flag) 1 else 0
      )
      hp_store[[tag]][[i]] <- list(
        depth = tune$depth,
        trees = tune$trees,
        innerK = innerK
      )
      
      st_outer <- prep_std(Xtr, Xte, std_flag)
      
      # IMPORTANT FIX: use ytr_bin from OUTER TRAIN ONLY
      df_tr <- data.frame(y = ytr_bin, as.data.frame(st_outer$Xtr))
      
      set.seed(seed + 7777 * i + if (std_flag) 1 else 0)
      
      fit <- gbm::gbm(
        y ~ .,
        data = df_tr,
        distribution = "bernoulli",
        n.trees = max(gbm_trees_grid),
        interaction.depth = tune$depth,
        shrinkage = gbm_shrinkage,
        n.minobsinnode = gbm_minobs,
        bag.fraction = gbm_bagfrac,
        train.fraction = 1.0,
        cv.folds = 0,
        verbose = FALSE
      )
      
      p_te <- predict(fit, newdata = as.data.frame(st_outer$Xva),
                      n.trees = tune$trees, type = "response")
      pred_store[[tag]][i] <- ifelse(p_te >= 0.5, "Brown", "White")
      
      rel_keep <- importance_gbm_rel(fit, n.trees = tune$trees, feat_names = feat_outer)
      v_full <- rep(0, length(feat_names))
      names(v_full) <- feat_names
      v_full[feat_outer] <- rel_keep
      IMP[[tag]][i, ] <- v_full
    }
  }
  
  # ----------------------------
  # fold-level results
  # ----------------------------
  truth_chr <- as.character(y_fac)
  fold_results <- data.frame(ID = ids, truth = truth_chr, stringsAsFactors = FALSE)
  for (nm in names(pred_store)) {
    fold_results[[paste0(nm, "_pred")]] <- pred_store[[nm]]
    fold_results[[paste0(nm, "_correct")]] <- as.integer(pred_store[[nm]] == truth_chr)
  }
  write.csv(fold_results, file.path(out_dir, "loocv_fold_results.csv"), row.names = FALSE)
  
  # ----------------------------
  # overall metrics
  # ----------------------------
  metric_row <- function(pred_vec) {
    pr <- factor(pred_vec, levels = levels(y_fac))
    met <- clf_metrics(y_fac, pr, positive_class = "Brown")
    c(
      acc = met$acc,
      bal_acc = met$bal_acc,
      sens = met$sensitivity,
      spec = met$specificity
    )
  }
  
  comparison <- do.call(rbind, lapply(names(pred_store), function(nm) {
    m <- metric_row(pred_store[[nm]])
    data.frame(
      model = nm,
      accuracy = unname(m["acc"]),
      balanced_accuracy = unname(m["bal_acc"]),
      sensitivity_brown = unname(m["sens"]),
      specificity_white = unname(m["spec"]),
      stringsAsFactors = FALSE
    )
  }))
  comparison <- comparison[order(-comparison$balanced_accuracy, -comparison$accuracy), ]
  write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)
  
  # ----------------------------
  # confusion matrices
  # ----------------------------
  conf_list <- lapply(names(pred_store), function(nm) {
    pr <- factor(pred_store[[nm]], levels = levels(y_fac))
    met <- clf_metrics(y_fac, pr, positive_class = "Brown")
    as.data.frame.matrix(met$conf)
  })
  names(conf_list) <- names(pred_store)
  saveRDS(conf_list, file.path(out_dir, "confusion_matrices_by_model.rds"))
  
  # ----------------------------
  # importance summaries
  # ----------------------------
  imp_summaries <- list(
    lasso_nostd = summarise_signed_importance(IMP$lasso_nostd, top_k = svm_top_k),
    lasso_std   = summarise_signed_importance(IMP$lasso_std,   top_k = svm_top_k),
    l2log_nostd = summarise_signed_importance(IMP$l2log_nostd, top_k = svm_top_k),
    l2log_std   = summarise_signed_importance(IMP$l2log_std,   top_k = svm_top_k),
    svm_lin_nostd = summarise_signed_importance(IMP$svm_lin_nostd, top_k = svm_top_k),
    svm_lin_std   = summarise_signed_importance(IMP$svm_lin_std,   top_k = svm_top_k),
    svm_poly_nostd = summarise_unsigned_importance(IMP$svm_poly_nostd, top_k = svm_top_k),
    svm_poly_std   = summarise_unsigned_importance(IMP$svm_poly_std,   top_k = svm_top_k),
    nb_nostd = summarise_unsigned_importance(IMP$nb_nostd, top_k = nb_top_k),
    nb_std   = summarise_unsigned_importance(IMP$nb_std,   top_k = nb_top_k),
    cart_nostd = summarise_unsigned_importance(IMP$cart_nostd, top_k = cart_top_k),
    cart_std   = summarise_unsigned_importance(IMP$cart_std,   top_k = cart_top_k),
    gbm_nostd = summarise_unsigned_importance(IMP$gbm_nostd, top_k = gbm_top_k),
    gbm_std   = summarise_unsigned_importance(IMP$gbm_std,   top_k = gbm_top_k)
  )
  
  safe_to_raw <- setNames(name_map$raw, name_map$safe)
  
  for (nm in names(imp_summaries)) {
    tmp <- imp_summaries[[nm]]
    if ("term" %in% names(tmp)) {
      tmp <- add_raw_names(tmp, safe_to_raw, "term")
      tmp$feature_label <- pretty_feature_name(tmp$raw_term)
    }
    write.csv(tmp, file.path(out_dir, paste0(nm, "_importance_summary.csv")), row.names = FALSE)
    imp_summaries[[nm]] <- tmp
  }
  
  # ----------------------------
  # feature-usage summary
  # ----------------------------
  usage_summary_df <- build_usage_summary(
    IMP = IMP, name_map = name_map, out_dir = out_dir,
    tau = 0.8, alpha = alpha_sig
  )
  
  # ----------------------------
  # number of "significant" features
  # ----------------------------
  signed_methods <- c("lasso_nostd","lasso_std","l2log_nostd","l2log_std","svm_lin_nostd","svm_lin_std")
  unsigned_methods <- c("svm_poly_nostd","svm_poly_std","nb_nostd","nb_std","cart_nostd","cart_std","gbm_nostd","gbm_std")
  
  sig_counts <- data.frame(
    model = c(signed_methods, unsigned_methods),
    n_sig = c(
      sapply(signed_methods, function(m) count_sig_features(IMP[[m]], signed = TRUE,  alpha = alpha_sig)),
      sapply(unsigned_methods, function(m) count_sig_features(IMP[[m]], signed = FALSE, alpha = alpha_sig))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(sig_counts, file.path(out_dir, "n_significant_features.csv"), row.names = FALSE)
  
  # ----------------------------
  # Save everything
  # ----------------------------
  saveRDS(
    list(
      seed = seed,
      fold_results = fold_results,
      comparison = comparison,
      confusion_matrices = conf_list,
      feature_name_map = name_map,
      hyperparameters_by_fold = hp_store,
      importance_matrices = IMP,
      importance_summaries = imp_summaries,
      usage_summary = usage_summary_df,
      n_significant_features = sig_counts,
      config = list(
        innerK_base = innerK_base,
        svm_cost_grid = svm_cost_grid,
        svm_poly_cost_grid = svm_poly_cost_grid,
        svm_poly_degree_grid = svm_poly_degree_grid,
        svm_poly_coef0_grid = svm_poly_coef0_grid,
        nb_prior_grid = nb_prior_grid,
        nb_laplace_grid = nb_laplace_grid,
        cart_cp_grid = cart_cp_grid,
        gbm_depth_grid = gbm_depth_grid,
        gbm_trees_grid = gbm_trees_grid,
        gbm_shrinkage = gbm_shrinkage,
        gbm_minobs = gbm_minobs,
        gbm_bagfrac = gbm_bagfrac,
        svm_poly_nperm = svm_poly_nperm,
        alpha_sig = alpha_sig
      )
    ),
    file.path(out_dir, "loocv_all_objects.rds")
  )
  
  list(
    comparison = comparison,
    sig_counts = sig_counts,
    out_dir = out_dir
  )
}

# =============================================================================
# MULTI-SEED RUNS
# =============================================================================
all_comp <- list()
all_sig  <- list()

for (s in seeds_to_run) {
  out_dir <- file.path(out_dir_base, paste0("seed_", s))
  res <- run_one_seed(seed = s, out_dir = out_dir, data_file = data_file)
  
  comp <- res$comparison
  comp$seed <- s
  all_comp[[as.character(s)]] <- comp
  
  sigc <- res$sig_counts
  sigc$seed <- s
  all_sig[[as.character(s)]] <- sigc
  
  cat("Done seed:", s, "->", out_dir, "\n")
}

comp_all <- dplyr::bind_rows(all_comp)
sig_all  <- dplyr::bind_rows(all_sig)

mean_metrics <- comp_all %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy   = sd(accuracy, na.rm = TRUE),
    mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    sd_balanced_accuracy   = sd(balanced_accuracy, na.rm = TRUE),
    mean_sensitivity_brown = mean(sensitivity_brown, na.rm = TRUE),
    sd_sensitivity_brown   = sd(sensitivity_brown, na.rm = TRUE),
    mean_specificity_white = mean(specificity_white, na.rm = TRUE),
    sd_specificity_white   = sd(specificity_white, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_balanced_accuracy), dplyr::desc(mean_accuracy))

sig_summary <- sig_all %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_n_sig = mean(n_sig, na.rm = TRUE),
    sd_n_sig   = sd(n_sig, na.rm = TRUE),
    min_n_sig  = min(n_sig, na.rm = TRUE),
    max_n_sig  = max(n_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_n_sig))

print(mean_metrics)
print(sig_summary)

write.csv(comp_all,     file.path(out_dir_base, "all_seeds_comparisons_long.csv"), row.names = FALSE)
write.csv(mean_metrics, file.path(out_dir_base, "mean_metrics_over_10_seeds.csv"), row.names = FALSE)
write.csv(sig_all,      file.path(out_dir_base, "all_seeds_n_significant_features_long.csv"), row.names = FALSE)
write.csv(sig_summary,  file.path(out_dir_base, "mean_n_significant_features_over_10_seeds.csv"), row.names = FALSE)

# =============================================================================
# GLOBAL FEATURE-DRIVER SUMMARIES ACROSS SEEDS
# =============================================================================
global_driver_objects <- build_global_feature_driver_summaries(
  out_dir_base = out_dir_base,
  seeds_to_run = seeds_to_run,
  top_k = global_top_k
)

# -----------------------------------------------------------------------------
# OPTIONAL: print top global drivers for key models
# -----------------------------------------------------------------------------
cat("\nTop global drivers for lasso_std:\n")
print(head(global_driver_objects$global_feature_summaries$lasso_std, 15))

cat("\nTop global drivers for svm_lin_std:\n")
print(head(global_driver_objects$global_feature_summaries$svm_lin_std, 15))

cat("\nTop consensus features across models:\n")
print(head(global_driver_objects$consensus_summary, 20))

cat("\nWrote:\n")
cat(" -", file.path(out_dir_base, "all_seeds_comparisons_long.csv"), "\n")
cat(" -", file.path(out_dir_base, "mean_metrics_over_10_seeds.csv"), "\n")
cat(" -", file.path(out_dir_base, "all_seeds_n_significant_features_long.csv"), "\n")
cat(" -", file.path(out_dir_base, "mean_n_significant_features_over_10_seeds.csv"), "\n")
cat(" -", file.path(out_dir_base, "GLOBAL_feature_consensus_long.csv"), "\n")
cat(" -", file.path(out_dir_base, "GLOBAL_feature_consensus_summary.csv"), "\n")
cat(" -", file.path(out_dir_base, "GLOBAL_feature_driver_objects.rds"), "\n")
cat(" - individual GLOBAL_feature_drivers_<model>.csv files in", out_dir_base, "\n")