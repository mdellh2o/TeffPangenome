#!/usr/bin/env bash
################
#Functional annotation#
################

for set in core_genes dispensable_genes unique_genes
do
#InterProScan
sh ./interproscan.sh -i ${set}.fasta -f tsv -goterms --cpu 28 -o interpro.${set}.txt # --> output: interpro.core_genes.txt, interpro.dispensable_genes.txt, interpro.unique_genes.txt
#Eggnog
emapper.py --data_dir ./eggnog-mapper-data -i ./${set}.fasta -o ${set}.eggnog --cpu 28 # --> output: core_genes.eggnog.emapper.annotations, dispensable_genes.eggnog.emapper.annotations, unique_genes.eggnog.emapper.annotations
done