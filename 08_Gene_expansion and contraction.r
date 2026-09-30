#----------------------------------------------------#
# GO on Quncho expanded and contracted gene families #
#----------------------------------------------------#

# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO_manual.html

setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/PostDoc/Teff/Teff_team/Manuscript/reviews/round_1/functional_annotation/GO_enrichment_exp_contr")
#BiocManager::install("topGO")
library(tidyverse)
library(topGO)
library(GO.db)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(patchwork)
library(viridis)


# Contracted and expanded gene 
Et_A_sig_contracted_genes <- read.delim("Et_A_sig_contracted_genes.tsv")
Et_B_sig_contracted_genes <- read.delim("Et_B_sig_contracted_genes.tsv")

Et_A_sig_expanded_genes <- read.delim("Et_A_sig_expanded_genes.tsv")
Et_B_sig_expanded_genes <- read.delim("Et_B_sig_expanded_genes.tsv")


# Quncho functional annotation - EggNog and InterProScan
quncho_egg <- read.delim("quncho.emapper.annotations", header=FALSE, comment.char="#")
quncho_egg <- setNames(data.frame(quncho_egg$V1, quncho_egg$V10), c("gene", "GO"))

quncho_ips <- read.delim("interpro_quncho.txt", header=FALSE)
quncho_ips <- setNames(data.frame(quncho_ips$V1, quncho_ips$V14), c("gene", "GO"))

##################
# 1. Merge functional annotation datasets 

# 1.1 Eggnog
quncho_egg <- data.frame(gene = quncho_egg$gene, GO = quncho_egg$GO)
quncho_egg$GO <- strsplit(as.character(quncho_egg$GO), ",") # list
quncho_egg <- quncho_egg[quncho_egg$GO != "-", ] 

# 1.2 InterProScan
quncho_ips <- data.frame(gene = quncho_ips$gene, GO = quncho_ips$GO)
quncho_ips <- quncho_ips %>%
  mutate(GO = str_replace_all(GO, "\\(.*?\\)", "")) %>% # remove (...)
  filter(GO != "-" & !is.na(GO)) %>% # filter out
  mutate(GO = strsplit(as.character(GO), "\\|")) %>% # list
  tidyr::unnest(GO) %>%
  mutate(GO = trimws(GO)) %>%      
  filter(GO != "" & !is.na(GO))  
quncho_ips <- quncho_ips %>%
  group_by(gene) %>%
  summarise(GO = paste(unique(GO), collapse = ","), .groups = "drop")
quncho_ips$GO <- strsplit(quncho_ips$GO, ",")

# 1.3 Merged Eggnog and InterProScan
all_genes <- full_join(quncho_egg, quncho_ips, by = "gene", suffix = c(".egg", ".ips"), relationship = "many-to-many")
all_genes <- all_genes %>%
  rowwise() %>%
  mutate(
    GO = list(unique(c(GO.egg, GO.ips)))
  ) %>%
  dplyr::select(gene, GO)


# 2. Extracted and contracted gene families
extract_genes <- function(df, group_name) {
  df %>%
    separate_rows(Gene_IDs, sep = ",\\s*") %>%
    transmute(gene = Gene_IDs, group = group_name)
}

exp_contr_genes <- rbind(
  extract_genes(Et_A_sig_contracted_genes, "SubgenomeA_contracted"),
  extract_genes(Et_B_sig_contracted_genes, "SubgenomeB_contracted"),
  extract_genes(Et_A_sig_expanded_genes,   "SubgenomeA_expanded"),
  extract_genes(Et_B_sig_expanded_genes,   "SubgenomeB_expanded")
)
exp_contr_genes <- exp_contr_genes[!is.na(exp_contr_genes[[1]]) & exp_contr_genes[[1]] != "(none)", ]
exp_contr_genes[[1]] <- sub("\\.1$", "", exp_contr_genes[[1]])


# 2.1 Merged All
all_genes <- left_join(all_genes, exp_contr_genes, by = "gene")



# 3. Create gene2GO list
gene2GO <- all_genes$GO
names(gene2GO) <- all_genes$gene

gene2GO <- lapply(gene2GO, unique) # remove duplicated genes
universe_genes <- names(gene2GO) # gene universe


# 4. Enrichment function
run_topgo_group <- function(target_group, ontology) {
  
  gene_vector <- ifelse(
    universe_genes %in% unique(all_genes$gene[all_genes$group == target_group]),1,0) # 1 for genes in targeted group
  names(gene_vector) <- universe_genes
  geneList <- factor(gene_vector)
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    annot = annFUN.gene2GO,
    nodeSize = 10, # to prune the GO hierarchy from the terms which have less than 10 annotated genes
    gene2GO = gene2GO)
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  table <- GenTable(
    GOdata,
    classicFisher = resultFisher,
    elimFisher = resultElim,
    orderBy = "elimFisher",
    topNodes = length(score(resultElim))
  )
  
  pvals <- score(resultElim)
  fdr_all <- p.adjust(pvals, method = "BH")
  table$FDR <- fdr_all[table$GO.ID]
  
  table <- table %>%
    filter(FDR < 0.05)
  
  if (nrow(table) == 0) {
    return(table)
  }
  
  table$Group <- target_group
  table$Ontology <- ontology
  table$FE <- table$Significant / table$Expected
  table$FDR_log10 <- -log10(table$FDR)
  
  return(table)
}


# 5. Run for the 3 groups
GO <- c("BP", "MF", "CC")
final_table <- data.frame()

for (ont in GO) {
  res_Et_A_contracted <- run_topgo_group(target_group = "SubgenomeA_contracted", ontology = ont)
  res_Et_B_contracted <- run_topgo_group(target_group = "SubgenomeB_contracted", ontology = ont)
  res_Et_A_expanded <- run_topgo_group(target_group = "SubgenomeA_expanded", ontology = ont)
  res_Et_B_expanded <- run_topgo_group(target_group = "SubgenomeB_expanded", ontology = ont)
  
  final_table <- rbind(final_table, res_Et_A_contracted, res_Et_B_contracted, res_Et_A_expanded, res_Et_B_expanded)
}

final_table$Term <- Term(GOTERM[final_table$GO.ID])

write.csv(final_table, "topGO_contr_exp.csv", row.names = FALSE)
#showSigOfNodes(res_core$GOdata, score(runTest(res_core$GOdata, algorithm = "elim", statistic = "fisher")), firstSigNodes = 5, useInfo = 'all')

final_table <- read.csv("topGO_contr_exp.csv")

# 6. plots
#groups <- c("Et_A_expanded", "Et_B_expanded") # "Et_A_contracted","Et_B_contracted"

#magma_dark <- viridis(256, option = "magma", direction = -1)[40:200]

final_table_p <- final_table %>%
  filter(FDR < 0.05) %>% # not sure if this is the best way to filter out !!!
  #mutate(Term = str_trunc(as.character(Term), width = 30)) %>%
  mutate(
    elimFisher = gsub("<", "", elimFisher),
    elimFisher = trimws(elimFisher),
    elimFisher = as.numeric(elimFisher)) 

df_ont <- final_table_p %>%
  filter(Ontology == "BP") # choose the ontology to plot 
# xlim_max <- max(df_ont$FE)  
# size_max <- max(df_ont$Significant)
#contracted
#groups <- c("Et_A_contracted", "Et_B_contracted") 
write.table(df_ont, 'GO_exp_contr.txt', quote=F, sep='\t', row.names = F)

groups <- c("SubgenomeA_contracted", "SubgenomeB_contracted") 
xlim_max <- c()  
size_max <- c()
global_max <- max(df_ont$Significant, na.rm = TRUE)

plot_list <- list()
heights <- c()
k <- 1

for (grp in groups) {#grp='SubgenomeB_contracted'
  
  df_sub <- df_ont %>%
    filter(Group == grp) %>%
    arrange(FDR) %>%
    slice_head(n = 10)
  
  show_legend <- (grp == "SubgenomeB_contracted")
  
  df_sub <- df_sub %>%
    mutate(Term_wrap = stringr::str_wrap(stringr::str_to_title(Term), width = 40))
    #mutate(Term_wrap = stringr::str_wrap(Term, width = 40))
  if (grp=='SubgenomeB_contracted'){
    df_sub$Term_wrap[2]<-c('Transcriptional Fidelity During\nTranscription Elongation By Rna Pol. II')
    df_sub$Term_wrap[3]<-c('Neg. Reg. Of Gene Expression\nVia Chromosomal Cpg Island Methylation')
  }
  xlim_max <- max(df_sub$FE)  
  size_max <-  max(df_sub$Significant)
  
  p <- ggplot(df_sub, aes(x = FE, y = fct_reorder(Term_wrap, FE), 
                          size = Significant, col = -log10(elimFisher))) +
    geom_point(alpha = 0.7) +
    #scale_color_gradientn(colors = magma_dark, name = "p-value") +
    scale_color_gradient(
      low = "#6BAED6",
      high =  "#082567", 
      name = "-log10(Pvalue)"
    ) +
    #scale_size_continuous(limits = c(1, size_max), range = c(1, 8), name = "Count") +
    #scale_size_continuous( range = c(1, 10), name = "Count") +
    scale_size_continuous(
      limits = c(0, global_max),
      breaks = c(10, 50, 100, 150),
      range = c(1, 10),
      name = "Count"
    ) +
    guides(size = guide_legend(override.aes = list(colour = "gray50"))) +
    coord_cartesian(xlim = c(1, xlim_max)) +
    theme_light() +
    ylab(grp) +
    xlab(ifelse(grp == "SubgenomeB_contracted", "Enrichment", "")) +
    theme(legend.position = if (show_legend) "right" else "none")
  
  
  plot_list[[k]] <- p
  heights[k] <- nrow(df_sub)
  k <- k + 1
}

wrap_plots(plot_list, ncol = 1, heights = heights)

pdf('GO_enrichment_contr_BP.pdf', width = 6 , height = 8)
wrap_plots(plot_list, ncol = 1, heights = heights)
dev.off()




#expanded
groups <- c("SubgenomeA_expanded", "SubgenomeB_expanded")
xlim_max <- c()  
size_max <- c()
global_max <- max(df_ont$Significant, na.rm = TRUE)

plot_list <- list()
heights <- c()
k <- 1

for (grp in groups) {#grp='SubgenomeB_expanded'
  
  df_sub <- df_ont %>%
    filter(Group == grp) %>%
    arrange(FDR) %>%
    slice_head(n = 10)
  
  show_legend <- (grp == "SubgenomeB_expanded")
  
  df_sub <- df_sub %>%
    mutate(Term_wrap = stringr::str_wrap(stringr::str_to_title(Term), width = 40))
  #mutate(Term_wrap = stringr::str_wrap(Term, width = 40))
  
  xlim_max <- max(df_sub$FE)  
  size_max <-  max(df_sub$Significant)
  
  p <- ggplot(df_sub, aes(x = FE, y = fct_reorder(Term_wrap, FE), 
                          size = Significant, col = -log10(elimFisher))) +
    geom_point(alpha = 0.7) +
    #scale_color_gradientn(colors = magma_dark, name = "p-value") +
    scale_color_gradient(
      low = "#6BAED6",
      high =  "#082567", 
      name = "-log10(Pvalue)"
    ) +
    #scale_size_continuous(limits = c(1, size_max), range = c(1, 8), name = "Count") +
    #scale_size_continuous( range = c(1, 10), name = "Count") +
    scale_size_continuous(
      limits = c(0, global_max),
      breaks = c(10, 50, 100, 150),
      range = c(1, 10),
      name = "Count"
    ) +
    guides(size = guide_legend(override.aes = list(colour = "gray50"))) +
    coord_cartesian(xlim = c(1, xlim_max)) +
    theme_light() +
    ylab(grp) +
    xlab(ifelse(grp == "SubgenomeB_expanded", "Enrichment", "")) +
    theme(legend.position = if (show_legend) "right" else "none")
  
  
  plot_list[[k]] <- p
  heights[k] <- nrow(df_sub)
  k <- k + 1
}

wrap_plots(plot_list, ncol = 1, heights = heights)

pdf('GO_enrichment_expanded_BP.pdf', width = 7 , height = 8)
wrap_plots(plot_list, ncol = 1, heights = heights)
dev.off()

getwd()

