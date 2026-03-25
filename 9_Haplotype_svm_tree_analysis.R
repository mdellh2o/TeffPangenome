# ============================================================
# Teff haplotype analysis
# ============================================================

# ----------------------------
# User parameters
# ----------------------------
seed <- 3
alpha_sig <- 0.05
adjust_method <- "BH"
svm_cost <- 1
svm_top_k <- 10

# Input files
heatmap_file <- "Quncho_4B_13150000_14016238.xlsx"
haplotype_file <- "haplo_onQuncho.txt"

# ----------------------------
# Packages
# ----------------------------
required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "viridis", "ggtext", "tibble",
  "data.table", "purrr", "e1071", "ggrepel", "scales", "grid"
)

missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggtext)
library(tibble)
library(data.table)
library(purrr)
library(e1071)
library(ggrepel)
library(scales)
library(grid)

set.seed(seed)

# ============================================================
# PART 1 - Chromosome 4B heatmap
# ============================================================

# ----------------------------
# 1) Phenotype mapping (Name -> phenotype)
# ----------------------------
name_vec <- c(
  "T224","T283","T365","T366","T87","Boni","T116","T330","T336","T-33",
  "T132","T288","T304","T379","T404","DZ","karadebi","Quncho","T177","T206",
  "T99","addisie","T412","T297","T345","Dabbi","Dtt2-02"
)

pheno_vec <- c(
  "White","White","White","White","White","White","Brown","White","Brown","White",
  "Brown","Brown","Brown","White","Brown","White","Brown","White","White","White",
  "White","White","Brown","Brown","Brown","Brown","White"
)

ph_map <- tibble(Name = name_vec, phenotype = pheno_vec)

# Colors for phenotype-specific y-axis labels
ph_cols <- c(
  "White" = "gold",
  "Brown" = "chocolate4"
)

# ----------------------------
# 2) Optional display-name changes for y labels
# ----------------------------
display_map <- c(
  "karadebi" = "Karadebi",
  "addisie"  = "Addisie",
  "T-33"     = "T33",
  "Dtt2-02"  = "Dtt2",
  "DZ"       = "DZ_01_354"
)

# ----------------------------
# 3) Read heatmap table
# ----------------------------
if (!file.exists(heatmap_file)) {
  stop("Heatmap file not found: ", heatmap_file)
}

tab <- read_excel(heatmap_file)

# ----------------------------
# 4) Explicit manual row order (top -> bottom)
# ----------------------------
manual_row_order_top_to_bottom <- c(
  "Quncho","T177","T379","DZ","addisie","T206","T365","T366","T412","T283",
  "T87","T330","T224","T99","T-33","Boni","T116","Dtt2-02","T304","T404",
  "T336","karadebi","Dabbi","T132","T288","T345","T297"
)

# Validate manual order against the Excel data
tab_names <- as.character(tab$Name)

if (anyDuplicated(manual_row_order_top_to_bottom)) {
  dupes <- unique(manual_row_order_top_to_bottom[duplicated(manual_row_order_top_to_bottom)])
  stop("Duplicate names in manual_row_order_top_to_bottom: ", paste(dupes, collapse = ", "))
}

missing_in_manual <- setdiff(tab_names, manual_row_order_top_to_bottom)
extra_in_manual   <- setdiff(manual_row_order_top_to_bottom, tab_names)

if (length(missing_in_manual) > 0) {
  stop(
    "These rows are in the Excel file but missing from manual_row_order_top_to_bottom: ",
    paste(missing_in_manual, collapse = ", ")
  )
}
if (length(extra_in_manual) > 0) {
  stop(
    "These names are in manual_row_order_top_to_bottom but NOT in the Excel file: ",
    paste(extra_in_manual, collapse = ", ")
  )
}

row_order_top_to_bottom <- manual_row_order_top_to_bottom

# ----------------------------
# 5) Reshape to long format and attach phenotype labels
# ----------------------------
tab_long <- tab %>%
  pivot_longer(
    cols = -Name,
    names_to = "position",
    values_to = "value"
  ) %>%
  mutate(
    position = as.character(position),
    value = as.numeric(value) / 20
  ) %>%
  left_join(ph_map, by = "Name") %>%
  mutate(
    Name = factor(Name, levels = rev(row_order_top_to_bottom))
  )

# ----------------------------
# 6) Build colored markdown labels for the y axis
# ----------------------------
ylabs_df <- tab_long %>%
  distinct(Name, phenotype) %>%
  mutate(
    Name_chr = as.character(Name),
    Name_display = recode(Name_chr, !!!display_map, .default = Name_chr),
    label_col = ifelse(is.na(phenotype), "black", unname(ph_cols[phenotype])),
    ylab = paste0("<span style='color:", label_col, ";'><b>", Name_display, "</b></span>")
  )

ylabs <- setNames(ylabs_df$ylab, ylabs_df$Name)

# ----------------------------
# 7) Helper functions for curly brackets on the right
# ----------------------------

# Convert a row span counted from the top into plot y coordinates
row_span_to_y <- function(from_top, to_top, n_rows, tile_half = 0.49) {
  stopifnot(from_top >= 1, to_top >= from_top, to_top <= n_rows)

  y_top_center <- n_rows - from_top + 1
  y_bot_center <- n_rows - to_top + 1

  c(
    y0 = y_bot_center - tile_half,  # lower edge
    y1 = y_top_center + tile_half   # upper edge
  )
}

# Build a smooth brace path that opens to the left
curly_brace_path <- function(x, y0, y1, width = 0.35, n = 250, id = "b1") {
  if (y1 < y0) {
    tmp <- y0
    y0 <- y1
    y1 <- tmp
  }

  t <- seq(0, 1, length.out = n)
  y <- y0 + (y1 - y0) * t
  x_offset <- (sin(2 * pi * t))^2
  x_path <- x - width * x_offset

  data.frame(
    x = x_path,
    y = y,
    brace_id = id,
    stringsAsFactors = FALSE
  )
}

# Build all braces from a specification table
build_braces <- function(spec, n_rows, x_base, default_width = 0.35, n = 250, y_gap = 0.12) {
  out <- vector("list", nrow(spec))

  for (i in seq_len(nrow(spec))) {
    yy <- row_span_to_y(spec$from_top[i], spec$to_top[i], n_rows)

    xo <- if ("x_offset" %in% names(spec)) spec$x_offset[i] else 0
    wd <- if ("width" %in% names(spec)) spec$width[i] else default_width
    col_i <- if ("color" %in% names(spec)) spec$color[i] else "grey20"

    # Apply the same vertical shrink to all brackets
    y0 <- yy["y0"] + y_gap
    y1 <- yy["y1"] - y_gap

    # Safety for one-row brackets
    if (y1 <= y0) {
      mid <- mean(c(yy["y0"], yy["y1"]))
      y0 <- mid - 0.05
      y1 <- mid + 0.05
    }

    tmp <- curly_brace_path(
      x = x_base + xo,
      y0 = y0,
      y1 = y1,
      width = wd,
      n = n,
      id = paste0("brace_", i)
    )
    tmp$brace_color <- col_i
    out[[i]] <- tmp
  }

  bind_rows(out)
}

# Build the text labels associated with each curly brace
build_brace_labels <- function(spec, n_rows, x_base, label_dx = 0.15, y_gap = 0.12) {
  if (!("label" %in% names(spec))) return(NULL)

  bind_rows(lapply(seq_len(nrow(spec)), function(i) {
    yy <- row_span_to_y(spec$from_top[i], spec$to_top[i], n_rows)
    xo <- if ("x_offset" %in% names(spec)) spec$x_offset[i] else 0

    txt_col <- if ("label_color" %in% names(spec)) {
      spec$label_color[i]
    } else if ("color" %in% names(spec)) {
      spec$color[i]
    } else {
      "grey20"
    }

    y0 <- yy["y0"] + y_gap
    y1 <- yy["y1"] - y_gap

    if (y1 <= y0) {
      mid <- mean(c(yy["y0"], yy["y1"]))
      y0 <- mid - 0.05
      y1 <- mid + 0.05
    }

    data.frame(
      x = x_base + xo + label_dx,
      y = mean(c(y0, y1)),
      label = spec$label[i],
      label_color = txt_col,
      stringsAsFactors = FALSE
    )
  }))
}

# ----------------------------
# 8) Numeric coordinates for custom plotting
# ----------------------------
pos_levels <- names(tab)[names(tab) != "Name"]
n_cols <- length(pos_levels)
n_rows <- length(row_order_top_to_bottom)

row_order_bottom_to_top <- rev(row_order_top_to_bottom)
y_breaks <- seq_len(n_rows)
y_labels <- unname(ylabs[row_order_bottom_to_top])

tab_long_num <- tab_long %>%
  mutate(
    x_id = match(position, pos_levels),
    y_id = n_rows - match(as.character(Name), row_order_top_to_bottom) + 1
  )

# ----------------------------
# 9) Curly bracket specification (rows counted from the top)
# ----------------------------
brace_spec <- tribble(
  ~from_top, ~to_top, ~x_offset, ~label,         ~color,          ~width, ~label_color,
  1,         11,      0.10,      "Haplotype 4",  "chartreuse3",   0.26,   "chartreuse3",
  12,        15,      0.10,      "Haplotype 1",  "#d95f02",       0.26,   "#d95f02",
  16,        18,      0.10,      "Haplotype 0",  "#1b9e77",       0.26,   "#1b9e77",
  19,        19,      0.10,      "Haplotype 6",  "azure4",        0.26,   "azure4",
  20,        20,      0.10,      "Haplotype 2",  "darkred",       0.26,   "darkred",
  21,        26,      0.10,      "Haplotype 3",  "darkgoldenrod", 0.26,   "darkgoldenrod",
  27,        27,      0.10,      "Haplotype 5",  "cyan3",         0.26,   "cyan3"
)

if (any(brace_spec$from_top < 1 | brace_spec$to_top > n_rows | brace_spec$from_top > brace_spec$to_top)) {
  stop("Invalid row spans in brace_spec. Rows must satisfy 1 <= from_top <= to_top <= n_rows.")
}

# ----------------------------
# 10) Build curly braces and labels
# ----------------------------
x_brace_base <- n_cols + 0.8
brace_gap <- 0.12

brace_paths <- build_braces(
  brace_spec,
  n_rows = n_rows,
  x_base = x_brace_base,
  default_width = 0.28,
  n = 250,
  y_gap = brace_gap
)

brace_labels <- build_brace_labels(
  brace_spec,
  n_rows = n_rows,
  x_base = x_brace_base,
  label_dx = 0.16,
  y_gap = brace_gap
)

# ----------------------------
# 11) Draw the heatmap
# ----------------------------
heatmap_plot <- ggplot(tab_long_num, aes(x = x_id, y = y_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.2, width = 0.98, height = 0.98) +
  scale_fill_viridis_c(option = "magma", name = "Distance") +
  scale_x_continuous(
    breaks = seq_len(n_cols),
    labels = pos_levels,
    expand = expansion(mult = c(0, 0.28))
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  geom_path(
    data = brace_paths,
    aes(x = x, y = y, group = brace_id, color = brace_color),
    inherit.aes = FALSE,
    linewidth = 0.9,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_text(
    data = brace_labels,
    aes(x = x, y = y, label = label, color = label_color),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  coord_cartesian(clip = "off") +
  labs(
    x = "Chr 4B sequence",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggtext::element_markdown(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

print(heatmap_plot)

# ============================================================
# PART 2 - Exploratory analysis + linear SVM
# ============================================================

# ----------------------------
# 12) Load the haplotype dataset
# ----------------------------
if (!file.exists(haplotype_file)) {
  stop("Haplotype data file not found: ", haplotype_file)
}

# Expected structure:
# - first column: sample ID
# - second column: phenotype
# - columns 3:9: haplotype markers
onquncho <- fread(haplotype_file, data.table = FALSE)

colnames(onquncho)[c(3:9)] <- c("Chr 1B","Chr 3B","Chr 4A","Chr 4B","Chr 5A","Chr 6A","Chr 9A")
rownames(onquncho) <- onquncho[, 1]
onquncho <- onquncho[, -1]

# Keep only the two phenotypes of interest and treat markers as factors
df <- onquncho %>%
  filter(phenotype %in% c("White", "Brown")) %>%
  mutate(phenotype = factor(phenotype, levels = c("White", "Brown"))) %>%
  mutate(across(-phenotype, ~ factor(as.character(.x))))

haplo_cols <- setdiff(names(df), "phenotype")

cat("\n============================\n")
cat("Dataset structure\n")
cat("============================\n")
str(df)

cat("\nPhenotype counts:\n")
print(table(df$phenotype))

# ----------------------------
# 13) Relevel each marker so that 0 is the reference, when present
# ----------------------------
set_ref0 <- function(d, haplo_cols) {
  d %>%
    mutate(across(
      all_of(haplo_cols),
      ~ if ("0" %in% levels(.x)) relevel(.x, ref = "0") else .x
    ))
}

df0 <- set_ref0(df, haplo_cols)

# ----------------------------
# 14) Plot raw haplotype counts by phenotype
# ----------------------------
plot_raw_haplotype_counts <- function(
  d,
  haplo_cols,
  title = "Raw haplotype level counts by phenotype",
  fill_colors = c("Brown" = "chocolate4", "White" = "gold"),
  breaks_by = 2,
  bar_width = 0.7,
  x_expand_add = 0.25
) {
  d_long <- d %>%
    select(phenotype, all_of(haplo_cols)) %>%
    pivot_longer(
      cols = all_of(haplo_cols),
      names_to = "marker",
      values_to = "level"
    ) %>%
    mutate(
      phenotype = factor(phenotype, levels = c("White", "Brown")),
      level_chr = as.character(level)
    )

  # Put numeric levels first, then non-numeric labels such as "Other"
  lvl_all <- unique(d_long$level_chr)
  num_all <- suppressWarnings(as.numeric(lvl_all))
  lvl_num <- lvl_all[!is.na(num_all)][order(num_all[!is.na(num_all)])]
  lvl_nonnum <- sort(lvl_all[is.na(num_all)])

  d_long <- d_long %>%
    mutate(level = factor(level_chr, levels = c(lvl_num, lvl_nonnum)))

  ymax <- d_long %>%
    count(marker, level) %>%
    summarise(mx = max(n), .groups = "drop") %>%
    pull(mx)

  ggplot(d_long, aes(x = level, fill = phenotype)) +
    geom_bar(width = bar_width) +
    facet_grid(. ~ marker, scales = "free_x", space = "free_x", switch = "x") +
    scale_x_discrete(expand = expansion(add = x_expand_add)) +
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = function(y) seq(0, ceiling(max(y)), by = breaks_by),
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "", y = "Haplotype counts", fill = "Phenotype", title = title) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", margin = margin(t = 4)),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.ticks.y = element_blank(),
      panel.spacing.x = unit(0.7, "lines"),
      plot.margin = margin(5.5, 5.5, 14, 5.5)
    )
}

count_plot <- plot_raw_haplotype_counts(
  df0,
  haplo_cols,
  title = "",
  bar_width = 0.65
)

print(count_plot)

# ----------------------------
# 15) Check associations among one-hot encoded levels
# ----------------------------
# We compute pairwise associations between dummy variables generated from
# the markers. Within-marker dummies are mutually exclusive by design, so
# we optionally keep only cross-marker comparisons.

term_to_marker <- function(term) sub("(Other|[0-9]+)$", "", term)

dummy_assoc_pairs_sig <- function(
  d,
  alpha = 0.05,
  adjust = "BH",
  cross_marker_only = TRUE
) {
  # One-hot encode all marker levels
  X <- model.matrix(phenotype ~ ., data = d)[, -1, drop = FALSE]
  terms <- colnames(X)

  cmb <- combn(seq_along(terms), 2)
  out <- vector("list", ncol(cmb))

  for (k in seq_len(ncol(cmb))) {
    i <- cmb[1, k]
    j <- cmb[2, k]

    xi <- X[, i]
    xj <- X[, j]

    tab_ij <- table(xi, xj)
    if (!all(dim(tab_ij) == c(2, 2))) {
      out[[k]] <- NULL
      next
    }

    # Correlation between two binary indicators is the phi coefficient
    corr_ij <- suppressWarnings(cor(xi, xj))
    p_ij <- suppressWarnings(stats::chisq.test(tab_ij, correct = FALSE)$p.value)

    out[[k]] <- data.frame(
      term1 = terms[i],
      term2 = terms[j],
      corr = as.numeric(corr_ij),
      p_value = p_ij,
      stringsAsFactors = FALSE
    )
  }

  res <- bind_rows(out) %>%
    mutate(
      abs_corr = abs(corr),
      p_adj = p.adjust(p_value, method = adjust),
      sig = p_adj < alpha
    )

  if (cross_marker_only) {
    res <- res %>%
      mutate(
        marker1 = term_to_marker(term1),
        marker2 = term_to_marker(term2)
      ) %>%
      filter(marker1 != marker2)
  }

  res %>% arrange(p_adj, desc(abs_corr))
}

pairsNC <- dummy_assoc_pairs_sig(
  df0,
  alpha = alpha_sig,
  adjust = adjust_method,
  cross_marker_only = TRUE
)

cat("\n============================\n")
cat("Top cross-marker dummy associations\n")
cat("============================\n")
print(head(pairsNC, 5))

# ============================================================
# PART 3 - Linear SVM with LOOCV and weight stability
# ============================================================

# ----------------------------
# 16) Helper: classification metrics
# ----------------------------
clf_metrics <- function(y_true, y_pred, positive_class = "Brown") {
  y_true <- factor(y_true)
  y_pred <- factor(y_pred, levels = levels(y_true))
  neg_class <- setdiff(levels(y_true), positive_class)[1]

  TP <- sum(y_true == positive_class & y_pred == positive_class)
  TN <- sum(y_true == neg_class     & y_pred == neg_class)
  FP <- sum(y_true == neg_class     & y_pred == positive_class)
  FN <- sum(y_true == positive_class & y_pred == neg_class)

  acc <- (TP + TN) / (TP + TN + FP + FN)
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)

  list(
    acc = acc,
    bal_acc = bal_acc,
    sensitivity = sens,
    specificity = spec,
    conf = table(truth = y_true, pred = y_pred)
  )
}

# ----------------------------
# 17) LOOCV linear SVM with fold-wise aligned weights
# ----------------------------
# Important idea:
# - each fold is standardized using the training set only
# - the linear SVM weight vector is extracted for each fold
# - the sign of the weights is aligned so that positive values point
#   toward the "Brown" class
# - weights are converted back to the original feature scale
svm_loocv_linear_aligned <- function(
  X,
  y,
  positive_class = "Brown",
  cost = 1,
  top_k = 10
) {
  X <- as.matrix(X)
  y <- droplevels(factor(y))

  stopifnot(positive_class %in% levels(y), nlevels(y) == 2)

  n <- nrow(X)
  p <- ncol(X)

  pred <- rep(NA_character_, n)
  dv_test <- rep(NA_real_, n)
  W <- matrix(0, nrow = n, ncol = p)
  colnames(W) <- colnames(X)

  topk_mat <- matrix(FALSE, nrow = n, ncol = p)
  colnames(topk_mat) <- colnames(X)

  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)

    Xtr <- X[tr, , drop = FALSE]
    Xte <- X[i,  , drop = FALSE]
    ytr <- y[tr]

    # Standardize within each training fold only
    sd_tr <- apply(Xtr, 2, sd)
    keep <- sd_tr > 0

    mu <- colMeans(Xtr[, keep, drop = FALSE])
    sigma <- sd_tr[keep]

    Xtr_s <- sweep(sweep(Xtr[, keep, drop = FALSE], 2, mu, "-"), 2, sigma, "/")
    Xte_s <- sweep(sweep(Xte[, keep, drop = FALSE], 2, mu, "-"), 2, sigma, "/")

    svm_i <- svm(
      x = Xtr_s,
      y = ytr,
      kernel = "linear",
      scale = FALSE,
      cost = cost,
      decision.values = TRUE
    )

    # Test prediction
    pr_te <- predict(svm_i, Xte_s, decision.values = TRUE)
    pred[i] <- as.character(pr_te)
    dv_test[i] <- as.numeric(attr(pr_te, "decision.values"))

    # Weight vector in standardized space
    w_scaled <- drop(t(svm_i$coefs) %*% svm_i$SV)

    # Align sign so that positive values indicate the positive class
    pr_tr <- predict(svm_i, Xtr_s, decision.values = TRUE)
    dv_tr <- as.numeric(attr(pr_tr, "decision.values"))
    mean_pos <- mean(dv_tr[ytr == positive_class], na.rm = TRUE)
    mean_neg <- mean(dv_tr[ytr != positive_class], na.rm = TRUE)
    sgn <- ifelse((mean_pos - mean_neg) >= 0, 1, -1)

    # Convert weights back to the original (unstandardized) feature scale
    w_orig <- w_scaled / sigma
    w_orig <- sgn * w_orig

    W[i, keep] <- w_orig

    # Store which terms are in the top-k by absolute weight
    k <- min(top_k, length(w_orig))
    top_idx <- order(abs(w_orig), decreasing = TRUE)[seq_len(k)]
    top_terms <- names(w_orig)[top_idx]
    topk_mat[i, colnames(X) %in% top_terms] <- TRUE
  }

  y_pred <- factor(pred, levels = levels(y))
  met <- clf_metrics(y, y_pred, positive_class = positive_class)

  # Aggregate fold-wise importance
  w_mean <- colMeans(W)
  w_sd <- apply(W, 2, sd)
  topk_freq <- colMeans(topk_mat)

  mean_sign <- sign(w_mean)
  sign_cons <- sapply(seq_len(p), function(j) {
    if (mean_sign[j] == 0) return(NA_real_)
    mean(sign(W[, j]) == mean_sign[j])
  })

  imp <- data.frame(
    term = colnames(X),
    mean_w = as.numeric(w_mean),
    sd_w = as.numeric(w_sd),
    mean_abs_w = abs(as.numeric(w_mean)),
    sign_consistency = as.numeric(sign_cons),
    topk_freq = as.numeric(topk_freq),
    stability_score = abs(as.numeric(w_mean)) / (as.numeric(w_sd) + 1e-8),
    direction = ifelse(
      w_mean > 0,
      positive_class,
      setdiff(levels(y), positive_class)[1]
    ),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(stability_score), desc(mean_abs_w))

  list(
    pred = y_pred,
    decision_values = dv_test,
    W = W,
    metrics = met,
    importance = imp,
    topk_mat = topk_mat
  )
}

# ----------------------------
# 18) Plot SVM stability
# ----------------------------
plot_svm_stability <- function(imp, top_n = 10, title = "SVM term stability") {
  pretty_term <- function(term) {
    term <- gsub("`", "", term)
    lvl  <- sub("^.*?(Other|[0-9]+)$", "\\1", term)
    mrk  <- sub("(Other|[0-9]+)$", "", term)
    mrk  <- gsub("\\.", " ", mrk)
    mrk  <- gsub("\\s+", " ", trimws(mrk))
    paste0(mrk, " = ", lvl)
  }

  imp2 <- imp %>%
    mutate(
      direction = ifelse(mean_w > 0, "Positive (Brown)", "Negative (White)"),
      label_pretty = vapply(term, pretty_term, character(1)),
      topk_freq = pmin(pmax(topk_freq, 0), 1)
    )

  top_terms <- imp2 %>%
    slice_max(order_by = stability_score, n = top_n, with_ties = FALSE) %>%
    pull(term)

  imp2 <- imp2 %>%
    mutate(label_to_show = ifelse(term %in% top_terms, label_pretty, ""))

  ggplot(imp2, aes(x = mean_abs_w, y = sd_w)) +
    geom_point(aes(color = direction, size = stability_score, alpha = topk_freq)) +
    scale_color_manual(values = c("Positive (Brown)" = "brown", "Negative (White)" = "gold")) +
    scale_size_continuous(name = "Stability score") +
    scale_alpha_continuous(
      name = "Top-k freq",
      range = c(0.10, 1),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    geom_label_repel(
      aes(label = label_to_show),
      fill = "white",
      color = "black",
      label.size = 0.25,
      label.padding = unit(0.15, "lines"),
      box.padding = 0.8,
      point.padding = 0.8,
      point.size = 6,
      force = 5,
      max.iter = 1e6,
      max.time = 5,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.size = 0.4,
      seed = 1,
      show.legend = FALSE
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      x = "|mean weight|",
      y = "sd(weight)",
      color = "Direction"
    ) +
    guides(
      color = guide_legend(order = 1, direction = "vertical"),
      size  = guide_legend(order = 2, direction = "vertical"),
      alpha = guide_legend(order = 3, direction = "vertical")
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.direction = "vertical",
      legend.box.just = "left",
      legend.title.align = 0,
      legend.text.align = 0,
      plot.margin = margin(5.5, 20, 5.5, 5.5)
    )
}

# ----------------------------
# 19) Build the design matrix and run SVM
# ----------------------------
# model.matrix() is kept because it is a very useful generic way to
# convert factor predictors into the one-hot encoded matrix required by
# the linear SVM.
X_svm <- model.matrix(phenotype ~ ., data = df0)[, -1, drop = FALSE]
y_svm <- df0$phenotype

svm_res <- svm_loocv_linear_aligned(
  X = X_svm,
  y = y_svm,
  positive_class = "Brown",
  cost = svm_cost,
  top_k = svm_top_k
)

cat("\n============================\n")
cat("SVM LOOCV confusion matrix\n")
cat("============================\n")
print(svm_res$metrics$conf)

cat("\nSVM LOOCV metrics:\n")
print(svm_res$metrics[c("acc", "bal_acc", "sensitivity", "specificity")])

cat("\nTop SVM terms by stability:\n")
print(head(svm_res$importance, 20))

svm_stability_plot <- plot_svm_stability(
  svm_res$importance,
  top_n = 5,
  title = ""
)

print(svm_stability_plot)

# ----------------------------
# 20) Compact interpretation helper
# ----------------------------
# Positive mean weight  -> evidence toward Brown
# Negative mean weight  -> evidence toward White
# Large |mean weight|   -> stronger average contribution
# Small sd(weight)      -> more stable across LOOCV folds
# Large topk_freq       -> often among the most influential terms
#
# Because the sample is very small (n = 27), these results should be
# described as stable within-sample patterns rather than guaranteed
# out-of-sample generalization.







################################################################################
# EXACT CHECK OF MINIMUM NUMBER OF HAPLOTYPE-LEVEL INDICATORS (binary features)
# THAT PERFECTLY CLASSIFY WHITE/BROWN ON THE OBSERVED n=27 SAMPLES

# Assumes:
# - df0 exists (your no-collapse dataframe with phenotype + haplo_cols)
# - haplo_cols exists (marker columns)
################################################################################

library(dplyr)

# ------------------------------------------------------------------------------
# 1) Build binary haplotype-level indicator columns (0/1), one per "marker = level"
#    Example column names: "Chr 4B_eq_3" or "H4B_eq_3"
# ------------------------------------------------------------------------------
make_level_binaries <- function(d, haplo_cols, drop_rare = 0, prefix = "") {
  stopifnot("phenotype" %in% names(d))
  stopifnot(all(haplo_cols %in% names(d)))
  
  bins <- lapply(haplo_cols, function(m) {
    x <- factor(d[[m]])
    tab <- table(x)
    keep_lvls <- names(tab)[tab > drop_rare]   # drop_rare=0 => keep ALL levels
    
    # one column per level (0/1)
    mat <- sapply(keep_lvls, function(lv) as.integer(x == lv))
    
    # if only one level remains, sapply can simplify to vector -> force matrix
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = 1)
      colnames(mat) <- keep_lvls
    }
    
    colnames(mat) <- paste0(prefix, m, "_eq_", keep_lvls)
    as.data.frame(mat, check.names = FALSE)
  })
  
  bind_cols(
    phenotype = droplevels(factor(d$phenotype)),
    bind_cols(bins)
  )
}

# Build binary matrix from your df0 + haplo_cols
df_bin_all <- make_level_binaries(df0, haplo_cols, drop_rare = 0)

# ------------------------------------------------------------------------------
# 2) OPTIONAL: If you already fitted a perfect rpart tree (res_full), inspect what
#    "5" means in that tree (depth vs #split nodes vs #distinct features)
# ------------------------------------------------------------------------------
inspect_rpart_usage <- function(tr) {
  fr <- tr$frame
  split_vars <- as.character(fr$var[fr$var != "<leaf>"])
  cat("Total split nodes:", length(split_vars), "\n")
  cat("Distinct split variables used:", length(unique(split_vars)), "\n")
  if (length(split_vars) > 0) {
    cat("Split variable frequency:\n")
    print(sort(table(split_vars), decreasing = TRUE))
  }
}

# Example (uncomment if res_full exists)
# inspect_rpart_usage(res_full$tree)

# ------------------------------------------------------------------------------
# 3) Robust converter: supports binary features coded as 0/1 OR YES/NO OR logical
# ------------------------------------------------------------------------------
prepare_Xy_from_df_bin <- function(df_bin, positive_class = "Brown") {
  stopifnot("phenotype" %in% names(df_bin))
  
  feat_names <- setdiff(names(df_bin), "phenotype")
  X_df <- df_bin[, feat_names, drop = FALSE]
  
  to01 <- function(z) {
    if (is.factor(z) || is.character(z)) {
      zz <- as.character(z)
      u <- unique(zz[!is.na(zz)])
      
      if (all(u %in% c("YES", "NO"))) return(as.integer(zz == "YES"))
      if (all(u %in% c("1", "0")))   return(as.integer(zz == "1"))
      
      stop("Unsupported character/factor coding in feature column. Values include: ",
           paste(head(u, 10), collapse = ", "))
    }
    
    if (is.logical(z)) return(as.integer(z))
    
    if (is.numeric(z) || is.integer(z)) {
      u <- unique(z[!is.na(z)])
      if (!all(u %in% c(0, 1))) {
        stop("Numeric feature is not binary 0/1. Values include: ",
             paste(head(u, 10), collapse = ", "))
      }
      return(as.integer(z))
    }
    
    stop("Unsupported feature type: ", paste(class(z), collapse = "/"))
  }
  
  X <- sapply(X_df, to01)
  X <- as.matrix(X)
  mode(X) <- "integer"
  colnames(X) <- feat_names
  
  y <- droplevels(factor(df_bin$phenotype))
  stopifnot(nlevels(y) == 2)
  if (!(positive_class %in% levels(y))) {
    stop("positive_class not found in phenotype levels.")
  }
  y01 <- as.integer(y == positive_class)
  
  list(X = X, y = y, y01 = y01, feat_names = feat_names)
}

xy <- prepare_Xy_from_df_bin(df_bin_all, positive_class = "Brown")
X <- xy$X
y <- xy$y
y01 <- xy$y01

cat("n samples =", nrow(X), "\n")
cat("p binary haplotype-level features =", ncol(X), "\n")
cat("Class counts:\n")
print(table(y))
cat("Sanity check (feature sums):\n")
print(summary(colSums(X)))

# ------------------------------------------------------------------------------
# 4) Exact test for a fixed subset of columns:
#    PERFECT iff every observed selected-feature pattern maps to one phenotype only
# ------------------------------------------------------------------------------
is_perfect_subset_idx <- function(cols, X, y01) {
  k <- length(cols)
  if (k == 0) return(length(unique(y01)) == 1L)
  
  # Bit-encode the k selected binary features into a pattern key (safe for k <= ~20)
  weights <- 2^(seq_len(k) - 1L)
  key <- as.integer(X[, cols, drop = FALSE] %*% weights)
  
  # Perfect if no key contains both classes
  groups <- split(y01, key)
  all(vapply(groups, function(v) length(unique(v)) == 1L, logical(1)))
}

# ------------------------------------------------------------------------------
# 5) Count ALL perfect subsets of size k (exact), optionally store the feature names
#    Uses combn(..., FUN=...) so it does NOT build the giant combinations matrix.
# ------------------------------------------------------------------------------
count_perfect_subsets_k <- function(X, y01, k,
                                    feat_names = colnames(X),
                                    keep_subsets = TRUE,
                                    max_store = Inf,
                                    verbose = TRUE) {
  p <- ncol(X)
  stopifnot(k >= 1, k <= p)
  
  total_combinations <- choose(p, k)
  if (verbose) {
    cat("\nCounting exact perfect subsets for k =", k, "\n")
    cat("Total combinations = choose(", p, ",", k, ") = ", format(total_combinations, scientific = FALSE), "\n", sep = "")
  }
  
  n_checked <- 0L
  n_perfect <- 0L
  hit_indices <- integer(0)       # positions in combn order
  hit_subsets_idx <- list()       # actual column index subsets
  hit_subsets_names <- list()     # feature names
  
  # callback called for each combination
  .fun <- function(cols) {
    n_checked <<- n_checked + 1L
    
    ok <- is_perfect_subset_idx(cols, X, y01)
    if (ok) {
      n_perfect <<- n_perfect + 1L
      
      if (keep_subsets && length(hit_subsets_idx) < max_store) {
        hit_indices <<- c(hit_indices, n_checked)
        hit_subsets_idx[[length(hit_subsets_idx) + 1L]] <<- cols
        hit_subsets_names[[length(hit_subsets_names) + 1L]] <<- feat_names[cols]
      }
    }
    
    # return scalar to keep combn happy; result is discarded
    FALSE
  }
  
  invisible(combn(p, k, FUN = .fun, simplify = TRUE))
  
  out <- list(
    k = k,
    p = p,
    n_combinations = total_combinations,
    n_checked = n_checked,
    n_perfect = n_perfect,
    perfect_indices = hit_indices,
    perfect_subsets_idx = hit_subsets_idx,
    perfect_subsets_names = hit_subsets_names
  )
  
  if (keep_subsets && length(hit_subsets_names) > 0) {
    out$perfect_subsets_df <- bind_rows(
      lapply(seq_along(hit_subsets_names), function(i) {
        data.frame(
          subset_id = i,
          t(hit_subsets_names[[i]]),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
    )
    colnames(out$perfect_subsets_df) <- c("subset_id", paste0("feature_", seq_len(k)))
  } else {
    out$perfect_subsets_df <- NULL
  }
  
  if (verbose) {
    cat("Checked combinations:", format(n_checked, scientific = FALSE), "\n")
    cat("Perfect combinations:", format(n_perfect, scientific = FALSE), "\n")
  }
  
  out
}

# ------------------------------------------------------------------------------
# 6) Find the MINIMUM k with at least one perfect subset (exact proof for k=1..k_max)
# ------------------------------------------------------------------------------
find_min_k_perfect <- function(X, y01,
                               feat_names = colnames(X),
                               k_max = 8,
                               keep_first_hits = TRUE,
                               max_store_first = Inf,
                               verbose = TRUE) {
  stopifnot(k_max >= 1)
  
  all_checked <- vector("list", k_max)
  
  for (k in seq_len(k_max)) {
    res_k <- count_perfect_subsets_k(
      X = X, y01 = y01, k = k,
      feat_names = feat_names,
      keep_subsets = keep_first_hits,
      max_store = max_store_first,
      verbose = verbose
    )
    
    all_checked[[k]] <- res_k
    
    if (res_k$n_perfect > 0) {
      return(list(
        min_k = k,
        first_nonzero = res_k,
        all_checked = all_checked[seq_len(k)]
      ))
    }
  }
  
  list(min_k = NA_integer_, first_nonzero = NULL, all_checked = all_checked)
}

# ------------------------------------------------------------------------------
# 7) Helper: inspect the pattern-to-class mapping induced by one perfect subset
#    (useful to understand the exact "criteria")
# ------------------------------------------------------------------------------
inspect_perfect_subset <- function(subset_features, X, y, feat_names = colnames(X)) {
  cols <- match(subset_features, feat_names)
  if (anyNA(cols)) stop("Some subset_features not found in feat_names.")
  
  Xs <- X[, cols, drop = FALSE]
  dfp <- as.data.frame(Xs)
  names(dfp) <- subset_features
  dfp$phenotype <- y
  
  # Pattern table
  pattern_tbl <- dfp %>%
    group_by(across(all_of(subset_features))) %>%
    summarise(
      n = n(),
      phenotypes = paste(sort(unique(as.character(phenotype))), collapse = "/"),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  
  pattern_tbl
}

# ------------------------------------------------------------------------------
# 8) RUN THE EXACT CHECKS
# ------------------------------------------------------------------------------

# (A) Confirm minimum k (exact)
# NOTE: k=7 can be computationally heavy. Start with k_max=6, then increase if needed.
min_res <- find_min_k_perfect(
  X = X, y01 = y01,
  feat_names = colnames(X),
  k_max = 8,              # you can lower to 6 first for speed
  keep_first_hits = TRUE,
  max_store_first = 1000, # cap how many perfect subsets to store
  verbose = TRUE
)

cat("\n========================================================\n")
cat("EXACT RESULT: minimum k =", min_res$min_k, "\n")
if (!is.null(min_res$first_nonzero)) {
  cat("Number of perfect combinations at k =", min_res$min_k, "is",
      min_res$first_nonzero$n_perfect, "\n")
}
cat("========================================================\n")

if (!is.null(min_res$first_nonzero) && !is.null(min_res$first_nonzero$perfect_subsets_df)) {
  cat("\nFirst 10 minimal perfect combinations:\n")
  print(head(min_res$first_nonzero$perfect_subsets_df, 10))
}

# (B) Specifically count ALL perfect combinations of size 5 (exact)
# This answers your direct question if k=5 is suspected.
res_k5 <- count_perfect_subsets_k(
  X = X, y01 = y01, k = 5,
  feat_names = colnames(X),
  keep_subsets = TRUE,
  max_store = Inf,     # set lower (e.g., 5000) if many hits and memory is a concern
  verbose = TRUE
)

cat("\nFor k = 5:\n")
cat("  total combinations =", format(res_k5$n_combinations, scientific = FALSE), "\n")
cat("  perfect combinations =", format(res_k5$n_perfect, scientific = FALSE), "\n")

if (!is.null(res_k5$perfect_subsets_df) && nrow(res_k5$perfect_subsets_df) > 0) {
  cat("\nFirst 20 perfect 5-feature combinations:\n")
  print(head(res_k5$perfect_subsets_df, 20))
}

# Optional save exact list
# if (!is.null(res_k5$perfect_subsets_df)) {
#   write.csv(res_k5$perfect_subsets_df, "perfect_5_feature_combinations.csv", row.names = FALSE)
# }

# ------------------------------------------------------------------------------
# 9) OPTIONAL: Compare with your tree result
#    If k=5 subsets are zero, that DOES NOT contradict a perfect depth-5 tree.
#    It means depth-5 tree separability exists, but no fixed set of 5 binary level
#    indicators alone is sufficient globally.
# ------------------------------------------------------------------------------

# Example: inspect the first perfect subset (if any)
if (!is.null(res_k5$perfect_subsets_names) && length(res_k5$perfect_subsets_names) > 0) {
  cat("\nPattern table for the FIRST perfect 5-feature subset:\n")
  print(inspect_perfect_subset(res_k5$perfect_subsets_names[[1]], X, y))
}

################################################################################
# END
################################################################################





library(dplyr)
library(rpart)
library(rpart.plot)

# ------------------------------------------------------------------------------
# Helper: extract leaf rules from an rpart object (if you don't already have it)
# ------------------------------------------------------------------------------
extract_leaf_rules <- function(tr) {
  fr <- tr$frame
  rn <- rownames(fr)
  
  leaf_nodes_chr <- rn[fr$var == "<leaf>"]
  leaf_nodes_num <- suppressWarnings(as.numeric(leaf_nodes_chr))
  rules_list <- path.rpart(tr, nodes = leaf_nodes_num, print.it = FALSE)
  
  out <- vector("list", length(leaf_nodes_chr))
  for (i in seq_along(leaf_nodes_chr)) {
    node_chr <- leaf_nodes_chr[i]
    ix <- match(node_chr, rn)
    
    leaf_id <- suppressWarnings(as.numeric(node_chr))
    n_leaf  <- as.integer(fr$n[ix])
    
    yv <- suppressWarnings(as.integer(fr$yval[ix]))
    pred_class <- if (!is.null(tr$ylevels) && !is.na(yv) && yv >= 1 && yv <= length(tr$ylevels)) {
      tr$ylevels[yv]
    } else {
      NA_character_
    }
    
    rule_lines <- rules_list[[i]]
    rule_txt <- if (is.null(rule_lines) || length(rule_lines) == 0) "" else paste(rule_lines, collapse = " & ")
    rule_txt <- sub("^root\\s*&\\s*", "", rule_txt)
    
    out[[i]] <- data.frame(
      leaf = leaf_id, n = n_leaf, pred = pred_class, rule = rule_txt,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(out)
}

# ------------------------------------------------------------------------------
# Helper: fit + plot a tree restricted to one perfect subset
# ------------------------------------------------------------------------------
plot_tree_for_perfect_subset <- function(df_bin_all,
                                         subset_features,
                                         maxdepth_cap = 10,
                                         split = "information",
                                         main_title = NULL) {
  stopifnot("phenotype" %in% names(df_bin_all))
  stopifnot(all(subset_features %in% names(df_bin_all)))
  
  dsub <- df_bin_all[, c("phenotype", subset_features), drop = FALSE]
  dsub$phenotype <- droplevels(factor(dsub$phenotype))
  
  # Fit a fully-grown tree on only these features
  tr <- rpart(
    phenotype ~ .,
    data = dsub,
    method = "class",
    parms = list(split = split),
    control = rpart.control(
      cp = 0,
      minsplit = 2,
      minbucket = 1,
      maxdepth = maxdepth_cap,
      xval = 0
    )
  )
  
  pred <- predict(tr, type = "class")
  acc <- mean(pred == dsub$phenotype)
  
  cat("\nSelected features:\n")
  print(subset_features)
  cat("\nTree accuracy on full dataset =", acc, "\n")
  cat("Number of split nodes =", sum(tr$frame$var != "<leaf>"), "\n")
  cat("Distinct split vars used =", length(unique(as.character(tr$frame$var[tr$frame$var != '<leaf>']))), "\n")
  
  # Plot
  rpart.plot(
    tr,
    type = 2,
    extra = 104,           # class + probabilities + n
    under = TRUE,
    fallen.leaves = TRUE,
    faclen = 0,            # don't shorten factor labels
    cex = 0.8,
    main = ifelse(is.null(main_title),
                  paste0("rpart on perfect subset (k = ", length(subset_features), "), acc = ", round(acc, 3)),
                  main_title)
  )
  
  cat("\nLeaf rules:\n")
  print(extract_leaf_rules(tr))
  
  invisible(list(tree = tr, data = dsub, acc = acc))
}

# ------------------------------------------------------------------------------
# 1) See the 5 perfect k=4 subsets
#    (You already have this in min_res$first_nonzero$perfect_subsets_df)
# ------------------------------------------------------------------------------
print(min_res$first_nonzero$perfect_subsets_df)

# ------------------------------------------------------------------------------
# 2) Pick one subset and plot its tree
#    Option A: from the names list (simplest)
# ------------------------------------------------------------------------------
subset4_1 <- min_res$first_nonzero$perfect_subsets_names[[1]]

res_plot_1 <- plot_tree_for_perfect_subset(
  df_bin_all = df_bin_all,
  subset_features = subset4_1,
  maxdepth_cap = 10,
  split = "information",
  main_title = "Perfect tree using minimal subset #1 (k=4)"
)

# ------------------------------------------------------------------------------
# 3) Try the other 4 minimal subsets (optional)
# ------------------------------------------------------------------------------
for (i in seq_along(min_res$first_nonzero$perfect_subsets_names)) {
  cat("\n===============================\n")
  cat("MINIMAL PERFECT SUBSET", i, "\n")
  cat("===============================\n")
  plot_tree_for_perfect_subset(
    df_bin_all = df_bin_all,
    subset_features = min_res$first_nonzero$perfect_subsets_names[[i]],
    maxdepth_cap = 10,
    split = "information",
    main_title = paste("Perfect tree for minimal subset", i, "(k=4)")
  )
}

# ------------------------------------------------------------------------------
# 4) (Optional) Plot using a subset selected by row from the data frame
# ------------------------------------------------------------------------------
get_subset_features_from_dfrow <- function(perfect_df, row_id) {
  row <- perfect_df[perfect_df$subset_id == row_id, , drop = FALSE]
  if (nrow(row) != 1) stop("subset_id not found or duplicated.")
  feats <- as.character(unlist(row[1, grepl("^feature_", names(row)), drop = FALSE]))
  feats[!is.na(feats)]
}

# Example:
# subset4_row3 <- get_subset_features_from_dfrow(min_res$first_nonzero$perfect_subsets_df, row_id = 3)
# plot_tree_for_perfect_subset(df_bin_all, subset4_row3, main_title = "Perfect tree for subset_id = 3")