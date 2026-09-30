#-------------------------------------#
# DNA methylation and Gene expression #
#-------------------------------------#

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(viridis)
library(ggpubr)
library(emmeans)
library(multcomp)
library(rstatix)
library(FSA)
library(multcompView)


## 1. LOAD DATA

feat <- read.delim("feature_full_table_cat.txt", header=TRUE)
length(unique(feat$gene))#51669

feat$subgenom <- sub(".*_([0-9]+)([AB])_.*", "\\2", feat$gene)
feat$feature_class <- sub(".*\\|(exon|intron)_.*", "\\1", feat$feature)
feat$TE_presence <- ifelse(feat$feature_TE_coverage > 0, "TE+", "TE-")
feat$feature_methylation <- ifelse(feat$feature_methylation == -1, 0, feat$feature_methylation)


## 2. TE PRESENCE AND METHYLATION

summary_df <- feat %>%
  filter(category %in% c("core", "dispensable")) %>%
  filter(feature_class %in% c("exon", "intron")) %>%
  group_by(category, TE_presence, feature_class) %>%
  summarise(
    mean_CpG = mean(feature_methylation, na.rm = TRUE),
    sd_CpG   = sd(feature_methylation, na.rm = TRUE),
    n        = sum(!is.na(feature_methylation)),
    se       = sd_CpG / sqrt(n),
    ci95     = qt(0.975, df = n - 1) * se,
    .groups  = "drop"
  )

df_test <- feat %>%
  filter(category %in% c("core", "dispensable")) %>%
  filter(feature_class %in% c("exon", "intron")) %>%
  mutate(
    TE_safe = recode(TE_presence, "TE+" = "TEpos", "TE-" = "TEneg"),
    group = interaction(category, TE_safe, feature_class, sep = "_")
  )

kw <- kruskal.test(feature_methylation ~ group, data = df_test)
print(kw)
dunn <- dunnTest(feature_methylation ~ group, data = df_test, method = "bh")
pvals <- setNames(dunn$res$P.adj, gsub(" - ", "-", dunn$res$Comparison))
letters_vec <- multcompLetters(pvals)$Letters

letters_df <- data.frame(
  group = names(letters_vec),
  letter = letters_vec,
  row.names = NULL) %>%
  separate(group, into = c("category", "TE_safe", "feature_class"), sep = "_") %>%
  mutate(
    TE_presence = recode(TE_safe, "TEpos" = "TE+", "TEneg" = "TE-"))

letters_pos <- summary_df %>%
  left_join(letters_df, by = c("category", "TE_presence", "feature_class")) %>%
  mutate(y = mean_CpG + ci95 + 2)  

ggplot(summary_df, aes(x = TE_presence, y = mean_CpG, fill = feature_class)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.8) +
  geom_errorbar(
    aes(ymin = mean_CpG - ci95, ymax = mean_CpG + ci95),
    position = position_dodge(0.9),
    width = 0.2) +
  geom_text(
    data = letters_pos,
    aes(x = TE_presence, y = y, label = letter, group = feature_class),
    position = position_dodge(0.9),
    size = 4) +
  facet_wrap(~category) +
  scale_fill_viridis_d(
    option = "magma",
    begin = 0.15,
    end = 0.85) +
  theme_bw() +
  labs(
    x = "TE presence",
    y = "Mean CpG methylation (%)",
    fill = "Feature class")


## 3. GENE EXPRESSION

gene_expr <- feat %>%
  distinct(gene, gene_expression, category) %>%
  mutate(log_expr = log10(gene_expression + 1))

pw <- pairwise_wilcox_test(gene_expr, log_expr ~ category,
                           p.adjust.method = "BH") # Pairwise Wilcoxon test
pval_vec <- setNames(pw$p.adj, paste(pw$group1, pw$group2, sep = "-"))

letters_out <- multcompLetters(pval_vec)$Letters
letters_df <- data.frame(category = names(letters_out), letter = letters_out)
label_pos <- gene_expr %>%
  group_by(category) %>%
  summarise(y = max(log_expr) * 1.05)

letters_df <- letters_df %>%
  left_join(label_pos, by = "category")

ggplot(gene_expr, aes(x = category, y = log_expr, fill = category)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  geom_text(
    data = letters_df,
    aes(x = category, y = y, label = letter),
    inherit.aes = FALSE,
    size = 5, fontface = "bold"
  ) +
  scale_fill_viridis_d(option = "magma", begin = 0.1, end = 0.9) +
  theme_bw() +
  labs(x = "Category", y = "log10(Expression + 1)")


## 4. METHYLATION vs EXPRESSION

feat_exon <- feat %>% filter(feature_class == "exon")

df_reg <- feat %>%
  distinct(gene, gene_expression, category) %>%
  left_join(
    feat_exon %>%
      group_by(gene) %>%
      summarise(
        mean_methylation = sum(feature_methylation * feature_length) / sum(feature_length),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  mutate(
    log_expr = log10(gene_expression + 1)
  )

df_plot <- df_reg %>%
  filter(category %in% c("core", "dispensable")) %>%
  group_by(category) %>%
  mutate(
    meth_group = ntile(mean_methylation, 4)
  ) %>%
  group_by(category, meth_group) %>%
  mutate(
    group_mean_meth = mean(mean_methylation, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    gene_class = factor(
      category,
      levels = c("core", "dispensable"),
      labels = c("Core genes", "Dispensable genes")
    ),
    group = interaction(meth_group, gene_class, sep = "_")
  )

kruskal.test(log_expr ~ group, data = df_plot)
dunn <- dunnTest(log_expr ~ group, data = df_plot, method = "bh")
pvals <- setNames(dunn$res$P.adj, gsub(" - ", "-", dunn$res$Comparison))

letters_vec <- multcompLetters(pvals)$Letters
letters_df <- data.frame(
  group = names(letters_vec),
  letter = letters_vec)
letters_pos <- df_plot %>%
  group_by(group, meth_group, gene_class) %>%
  summarise(
    y = max(log_expr) + 0.2,
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "group")

ggplot(df_plot,
       aes(x = factor(meth_group),
           y = log_expr,
           fill = gene_class)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.3) +
  geom_text(
    data = letters_pos,
    aes(x = factor(meth_group), y = y, label = letter, group = gene_class),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_viridis_d(
    option = "magma",
    begin = 0.1,
    end = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Methylation quartiles",
    y = "Expression (log10)",
    fill = "Gene class"
  )


####
# GROUPING ALL GENES

df_reg <- feat %>%
  distinct(gene, gene_expression, category) %>%
  left_join(
    feat %>%
      group_by(gene) %>%
      summarise(
        mean_methylation = sum(feature_methylation * feature_length) / sum(feature_length),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  mutate(
    log_expr = log10(gene_expression + 1)
  )

df_plot <- df_reg %>%
  filter(category %in% c("core", "dispensable")) %>%
  mutate(
    meth_group = ntile(mean_methylation, 4)) %>%
  group_by(meth_group) %>%
  mutate(
    group_mean_meth = mean(mean_methylation, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    group = factor(meth_group))

kruskal.test(log_expr ~ group, data = df_plot)
dunn <- dunnTest(log_expr ~ group, data = df_plot, method = "bh")
pvals <- setNames(dunn$res$P.adj, gsub(" - ", "-", dunn$res$Comparison))
letters_vec <- multcompLetters(pvals)$Letters
letters_df <- data.frame(
  group = names(letters_vec),
  letter = letters_vec)

letters_pos <- df_plot %>%
  group_by(group, meth_group) %>%
  summarise(
    y = max(log_expr) + 0.2,
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "group")

axis_labels <- df_plot %>%
  distinct(meth_group, group_mean_meth) %>%
  arrange(meth_group) %>%
  mutate(label = sprintf("%.2f", group_mean_meth))

ggplot(df_plot,
       aes(x = factor(meth_group),
           y = log_expr)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, fill = "darkviolet") +
  geom_text(
    data = letters_pos,
    aes(x = factor(meth_group), y = y, label = letter),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_x_discrete(
    breaks = axis_labels$meth_group,
    labels = axis_labels$label
  ) +
  theme_bw() +
  labs(
    x = "Methylation level (mean per quartile)",
    y = "Expression (log10)"
  )
