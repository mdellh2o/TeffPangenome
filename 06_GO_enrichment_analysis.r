#-----------------------#
#GO Enrichment Analysis#
#-----------------------#

#GO term enrichment of core, dispensable and unique genes, 
#based on the eggNOG-mapper and InterProScan annotations of the orthogroup

# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO_manual.html

#BiocManager::install("topGO")
# install.packages('tidyverse')
# install.packages('ggvenn')
library(tidyverse)
library(topGO)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(patchwork)
library(viridis)
library(Rgraphviz)
library(ggvenn)


# EggNog
setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/PostDoc/Teff/Teff_team/Manuscript/reviews/round_1/eggnog_out")
core_egg <- read.delim("single_core_in_og.eggnog.emapper.annotations", header=FALSE, comment.char="#")
core_egg <- setNames(data.frame(core_egg$V1, core_egg$V10), c("gene", "GO"))
dispensable_egg <- read.delim("single_disp_in_og.eggnog.emapper.annotations", header=FALSE, comment.char="#")
dispensable_egg <- setNames(data.frame(dispensable_egg$V1, dispensable_egg$V10), c("gene", "GO"))
unique_egg <- read.delim("single_uniq_in_og.eggnog.emapper.annotations", header=FALSE, comment.char="#")
unique_egg <- setNames(data.frame(unique_egg$V1, unique_egg$V10), c("gene", "GO"))


# InterProScan
setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/PostDoc/Teff/Teff_team/Manuscript/reviews/round_1/interpro_out")
core_ips <- read.delim("interpro.single_core_in_og.txt", header=FALSE)
core_ips <- setNames(data.frame(core_ips$V1, core_ips$V14), c("gene", "GO"))
dispensable_ips <- read.delim("interpro.single_disp_in_og.txt", header=FALSE)
dispensable_ips <- setNames(data.frame(dispensable_ips$V1, dispensable_ips$V14), c("gene", "GO"))
unique_ips <- read.delim("interpro.single_uniq_in_og.txt", header=FALSE)
unique_ips <- setNames(data.frame(unique_ips$V1, unique_ips$V14), c("gene", "GO"))

##################
# 1. Merge datasets 
# 1.1 Eggnog
all_genes_egg <- rbind(
  data.frame(gene = core_egg$gene, GO = core_egg$GO, group = "core"),
  data.frame(gene = dispensable_egg$gene, GO = dispensable_egg$GO, group = "dispensable"),
  data.frame(gene = unique_egg$gene, GO = unique_egg$GO, group = "unique"))

all_genes_egg$GO <- strsplit(as.character(all_genes_egg$GO), ",") # list
all_genes_egg <- all_genes_egg[all_genes_egg$GO != "-", ] 

# all_genes_egg <- all_genes_egg %>%
#   filter(GO != "-")
# 
# all_genes_egg$GO <- strsplit(as.character(all_genes_egg$GO), ",")


# 1.2 InterProScan
all_genes_ips <- rbind(
  data.frame(gene = core_ips$gene, GO = core_ips$GO, group = "core"),
  data.frame(gene = dispensable_ips$gene, GO = dispensable_ips$GO, group = "dispensable"),
  data.frame(gene = unique_ips$gene, GO = unique_ips$GO, group = "unique"))

all_genes_ips <- all_genes_ips %>%
  mutate(GO = str_replace_all(GO, "\\(.*?\\)", "")) %>% # remove (...)
  filter(GO != "-" & !is.na(GO)) %>% # filter out
  mutate(GO = strsplit(as.character(GO), "\\|")) %>% # list
  tidyr::unnest(GO) %>%
  mutate(GO = trimws(GO)) %>%      
  filter(GO != "" & !is.na(GO))    

all_genes_ips <- all_genes_ips %>%
  group_by(gene, group) %>%
  summarise(GO = paste(unique(GO), collapse = ","), .groups = "drop")
all_genes_ips$GO <- strsplit(all_genes_ips$GO, ",")


# # 2. chose which functional annotation program to use
# 
# all_genes <- all_genes_egg # only eggnog
# all_genes <- all_genes_ips # only interproscan

# merged
all_genes <- full_join(all_genes_egg, all_genes_ips, by = "gene", suffix = c(".egg", ".ips"), relationship = "many-to-many")
all_genes <- all_genes %>%
  rowwise() %>%
  filter(
    is.na(group.egg) | is.na(group.ips) | group.egg == group.ips
  ) %>%
  mutate(
    GO = list(unique(c(GO.egg, GO.ips))),
    group = ifelse(!is.na(group.egg), group.egg, group.ips)
  ) %>%
  ungroup() %>%
  dplyr::select(gene, GO, group)


# 3. Create gene2GO list
gene2GO <- all_genes$GO
names(gene2GO) <- all_genes$gene

gene2GO <- lapply(gene2GO, unique) # remove duplicated genes
universe_genes <- names(gene2GO) # gene universe


# 4. Enrichment function
run_topgo_group <- function(target_group, ontology) {
  
  #target_group='core'
  #ontology='BP'
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
  res_core <- run_topgo_group(target_group = "core", ontology = ont)
  res_disp <- run_topgo_group(target_group = "dispensable", ontology = ont)
  res_unique <- run_topgo_group(target_group = "unique", ontology = ont)
  
  final_table <- rbind(final_table, res_core, res_disp, res_unique)
}


#write.csv(final_table, "topGO_interproscan_top20.csv", row.names = FALSE)
#showSigOfNodes(res_core$GOdata, score(runTest(res_core$GOdata, algorithm = "elim", statistic = "fisher")), firstSigNodes = 5, useInfo = 'all')


# 6. plots
groups <- c("core","dispensable","unique")
magma_dark <- viridis(256, option = "magma", direction = -1)[40:200]

final_table_p <- final_table %>%
  filter(FDR < 0.05) %>%
  mutate(Term = str_trunc(as.character(Term), width = 30)) %>%
  mutate(
    elimFisher = gsub("<", "", elimFisher),
    elimFisher = trimws(elimFisher),
    elimFisher = as.numeric(elimFisher))
write.table(final_table_p, 'final_table_p.txt', quote=F, sep='\t', row.names = F)

#library(data.table)
#setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/PostDoc/Teff/Teff_team/Manuscript/reviews/round_1/functional_annotation/GO_enrichment_v2")
#final_table_p<-fread('final_table_p_long.txt', data.table=F)

df_ont <- final_table_p %>%
  filter(Ontology == "BP") # choose the ontology to plot 
# xlim_max <- max(df_ont$FE)  
# size_max <- max(df_ont$Significant)

plot_list <- list()
heights <- c()
k <- 1
xlim_max <- c()  
size_max <- c()
global_max <- max(df_ont$Significant, na.rm = TRUE)

library(stringr)

for (grp in groups) {#grp='unique'
  
  df_sub <- df_ont %>%
    filter(Group == grp)
  
  show_legend <- (grp == "dispensable")
  
  df_sub <- df_sub %>%
    mutate(Term_wrap = stringr::str_wrap(stringr::str_to_title(Term), width = 40)) 
  # %>%
  #   mutate(sig = -log10(elimFisher))
  
  xlim_max <- max(df_sub$FE)  
  size_max <-  max(df_sub$Significant)
  
  p <- ggplot(df_sub, aes(x = FE, y = fct_reorder(Term_wrap, FE), 
                          size = Significant, col = -log10(elimFisher))) +
    geom_point(alpha = 0.7) +
    #scale_color_gradientn(colors = magma_dark,  name = "p-value") +
    # scale_color_gradientn(
    #   colours = c(
    #     "#081D58", # very dark blue
    #     "#225EA8",
    #     "#41B6C4",
    #     "#C7E9B4",
    #     "#FFFFD9" # very light
    #     ),
    #   name = "p-value"
    # ) +
    scale_color_gradient(
      low = "#6BAED6",
      high =  "#082567", 
      name = "-log10(Pvalue)"
    ) +
    #scale_size_continuous(limits = c(1, size_max), range = c(1, 8), name = "Count") +
    #scale_size_continuous( range = c(1, 10), name = "Count") +
    scale_size_continuous(
      limits = c(5, global_max),
      breaks = c(5, 50, 500, 5000),
      range = c(1, 10),
      name = "Count"
    ) +
    # scale_size_area(
    #   max_size = 10,
    #   breaks = c(5, 50, 500, 5000),
    #   limits = c(0, global_max),
    #   name = "Count"
    # ) +
    guides(size = guide_legend(override.aes = list(colour = "gray50"))) +
    coord_cartesian(xlim = c(1, xlim_max)) +
    theme_light() +
    ylab(grp) +
    xlab(ifelse(grp == "unique", "Enrichment", "")) +
    theme(legend.position = if (show_legend) "right" else "none")
  
  write.table(df_sub, paste0(grp,'_dfBP.txt'), quote=F, sep='\t', row.names = F)
  
  plot_list[[k]] <- p
  heights[k] <- nrow(df_sub)
  k <- k + 1
}

wrap_plots(plot_list, ncol = 1, heights = heights)


pdf('GO_enrichment_BP.pdf', width = 6 , height = 8)
wrap_plots(plot_list, ncol = 1, heights = heights)
dev.off()

#which genes have GO:0010200

#response to chitin
TorF<-c()
for (i in 1:length(gene2GO)){#i=2
  TorF[i]<-'GO:0010200' %in% gene2GO[[i]]
}
chitin_genes<-names(gene2GO[which(TorF==TRUE)])

#response to symbiotic fungus
TorF<-c()
for (i in 1:length(gene2GO)){#i=2
  TorF[i]<-'GO:0009610' %in% gene2GO[[i]]
}
symfun_genes<-names(gene2GO[which(TorF==TRUE)])

#######
unique(symfun_genes,chitin_genes)
intersect(symfun_genes,chitin_genes)

sum(grepl("_[0-9]+A_", symfun_genes)) # number of A --> 22

sum(grepl("_[0-9]+B_", symfun_genes)) # number of B --> 23

# Extract the number between the first and second underscore
numA <- as.integer(sub("^[^_]+_([0-9]+)[A]_.*$", "\\1", symfun_genes))
numB <- as.integer(sub("^[^_]+_([0-9]+)[B]_.*$", "\\1", symfun_genes))
num <- as.integer(sub("^[^_]+_([0-9]+)[AB]_.*$", "\\1", symfun_genes))
# Count occurrences
table(numA)
table(numB)
table(num)


# num
# 1  2  3  4  6  7  8  9 
# 10 10  2  9  2  2  6  4 
