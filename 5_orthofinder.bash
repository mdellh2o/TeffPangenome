#!/usr/bin/env bash
##############
#Orthofinder #
##############


for file in ../Helixer_annotation/*.gff3
do
	ass=$(basename $file .filt.EDTA.gff3)

	ref=../${ass}*_chr.filt.fasta
	gff=${file}
	dir=../orthofinder/peptide/
	
	agat_sp_extract_sequences.pl \
	-g ${gff} \
	-f ${ref} \
	-p \
	-o ${dir}/${ass}.fa
	
	dir=../orthofinder/bed/
	awk 'OFS="\t" {split($NF,name,";");split(name[1],id,"="); if($3=="mRNA") print $1,$4,$5,id[2]}' ${gff} > ${dir}/${ass}.bed

done


##BUSCO analysis on filtered Helixer annotations


for fa in ../orthofinder/peptide/*.fa
do
	name=$(basename $fa .fa)

	busco \
	 -i ${fa} \
	 -m prot \
	 -l poales_odb12 \
	 -c 10 \
	 -o ${name} \
	 --out_path ../helixer_annotation/helixer_annotation_filt_EDTA_BUSCO \
	 --plot_percentages 

	busco \
	 --plot ../helixer_annotation/helixer_annotation_filt_EDTA_BUSCO/${name} 

done


## Run OrthoFinder (OrthoFinder version 2.5.5)

orthofinder \
 -f ../orthofinder/peptide/ \
 -X \
 -t 20

# Sequence search program [Default = diamond] 
# Merging single and multiple Orthogroups and Removing 1 and 2 orthogroups 
head -n 1 Orthogroups.tsv > orthogroups_all.tsv
tail -n +2 Orthogroups.tsv >> orthogroups_all.tsv
tail -n +2 Orthogroups_UnassignedGenes.tsv >> orthogroups_all.tsv


#############################################################################################
#####################################
#R script to split subgenome A and B#
#####################################


## Splitting A and B
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(readr)

orthogroups_all <- read.delim("../orthofinder/orthogroups_all.tsv", header=TRUE)
orthogroups_AB_all <- orthogroups_all %>% select(Orthogroup)

for (col in colnames(orthogroups_all)[-1]) {
  split_list <- str_split(orthogroups_all[[col]], ",\\s*")

  base_name <- str_extract(col, "^[^_]+")
  
  A <- sapply(split_list, function(x) paste(x[str_detect(x, "[0-9]+A_")], collapse=", "))
  B <- sapply(split_list, function(x) paste(x[str_detect(x, "[0-9]+B_")], collapse=", "))
 
  orthogroups_AB_all[[paste0(base_name, "_A")]] <- A
  orthogroups_AB_all[[paste0(base_name, "_B")]] <- B
}

write.table(orthogroups_AB_all,
            file = "../orthofinder/orthogroups_AB_all.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            na = "")

orthogroups_A_all <- orthogroups_AB_all[, c("Orthogroup", grep("_A$", colnames(orthogroups_AB_all), value=TRUE))]
orthogroups_B_all <- orthogroups_AB_all[, c("Orthogroup", grep("_B$", colnames(orthogroups_AB_all), value=TRUE))]


dfs <- list(orthogroups_A_all, orthogroups_B_all) # orthogroups_all
suffixes <- c("_A", "_B")

for (i in 1:2){
  
  df <- dfs[[i]]
  suffix <- suffixes[i]
  df <- df %>% mutate_all(~na_if(., ""))
  empty_rows <- apply(df[,-1], 1, function(x) all(is.na(x))) # removing empty rows
  df <- df[!empty_rows, ]

  # deviding genes for core, dispensable and unique
  core <- c()
  dispensable <- c()
  unique <- c() 

  for (j in 1:nrow(df)){
    empty <- sum(is.na(df[j, ]))
    if (empty == 0){
      core <- c(core, df[j,1])
    }
    else if (empty == (ncol(df)-2)){
      unique <- c(unique, df[j,1])
    }
    else{
      dispensable <- c(dispensable, df[j,1])
    }
  }
  write.table(core, paste0("../orthofinder/core_genes", suffix, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(dispensable, paste0("../orthofinder/dispensable_genes", suffix, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique, paste0("../orthofinder/unique_genes", suffix, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  
  ## Summary table
  summary_table <- data.frame(
    Sample = colnames(df)[-1],  
    Core = 0,
    Dispensable = 0,
    Unique = 0
  )

  for (col in 2:ncol(df)) {
    sample_name <- colnames(df)[col]

    sample_col <- df[, col]

    summary_table[summary_table$Sample == sample_name, "Core"] <-
      sum(df[,1] %in% core & !is.na(sample_col))

    summary_table[summary_table$Sample == sample_name, "Dispensable"] <-
      sum(df[,1] %in% dispensable & !is.na(sample_col))

    summary_table[summary_table$Sample == sample_name, "Unique"] <-
      sum(df[,1] %in% unique & !is.na(sample_col))
  }
  write.table(summary_table, paste0("../orthofinder/summary_table", suffix, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

#############################################################################################

##############
#back to bash#
##############
## Extracting core/dispensable/unique gene sequences

category="core" # core / dispensable / unique

sort -u \
 ../orthofinder/${category}_genes_A.txt \
 ../orthofinder/${category}_genes_B.txt \
| while read OG; do
    fasta="../Orthogroup_Sequences/${OG}.fa"
    [[ ! -f "$fasta" ]] && continue

    # If Quncho is present in the file, take only the first sequence of Quncho
    if grep -q "^>Quncho_" "$fasta"; then
        awk '/^>Quncho_/{print; flag=1; next} /^>/ && flag {exit} flag' "$fasta"
    else
        # Else take the first sequence of the file
        awk '/^>/{if(seen++) exit}1' "$fasta"
    fi
done \
 > ../orthofinder/${category}_genes.fasta

grep -c "^>" ../orthofinder/${category}_genes.fasta
