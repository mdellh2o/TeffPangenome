#!/usr/bin/env bash
## Local Pangenome Graph of the TT2 Locus (Chr 4B)
###############################################################################
# STEP 1: Perform all-vs-all genome alignment with wfmash
# Input: pangenome FASTA containing all assemblies
# Output: PAF file containing pairwise alignments between genomes
###############################################################################

conda activate pggb

cd /home/shared/shared_all/e.teff/genomes/teff/assembly/chr_genome

# PanSN-formatted all-vs-all alignment.
# -Y '#' tells wfmash to treat accessions sharing the same genome name
# as belonging to one assembly and exclude self-genome mappings.
wfmash -Y '#' E_tef_pangenome.fasta.gz -t 18 \
    > ./wfmash/E_tef_pangenome.paf

conda deactivate


###############################################################################
# STEP 2: Query homologous regions around the locus of interest with IMPG
# Reference locus:
# Quncho_chr#1#4B:13910000-13950000
#
# Output:
# BED file with corresponding intervals across all accessions
###############################################################################

conda activate impg

cd /home/shared/shared_all/e.teff/impg

impg query \
    -a /home/shared/shared_all/e.teff/genomes/teff/assembly/chr_genome/wfmash/E_tef_pangenome.paf \
    -r Quncho_chr#1#4B:13910000-13950000 \
    -d 100 \
    -x \
    -o bed \
    > tef_TT2_chr4B.bed

conda deactivate


###############################################################################
# STEP 3: Merge overlapping/nearby homologous intervals
#
# Sort BED intervals
# Merge intervals separated by ≤1000 bp
# Retain only regions larger than 10 kb
###############################################################################

conda activate bedtools

bedtools sort \
    -i tef_TT2_chr4B.bed \
    > tef_TT2_chr4B_sorted.bed

bedtools merge \
    -d 1000 \
    -i tef_TT2_chr4B_sorted.bed \
    > tef_TT2_chr4B_merged.bed

# Save merged intervals >10 kb in samtools region format
awk '$3 - $2 > 10000 { print $1":"$2"-"$3 }' \
    tef_TT2_chr4B_merged.bed \
    > tef_TT2_chr4B_merged.regions

conda deactivate


###############################################################################
# STEP 4: Extract interval
#
# Output:
# Compressed FASTA containing only the selected homologous regions
###############################################################################

conda activate samtools

samtools faidx \
    -r tef_TT2_chr4B_merged.regions \
    /home/shared/shared_all/e.teff/genomes/teff/assembly/chr_genome/E_tef_pangenome.fasta.gz \
    | bgzip -@24 -l9 \
    > tef_TT2_chr4B_merged.fa.gz

# Index extracted FASTA
samtools faidx tef_TT2_chr4B_merged.fa.gz

conda deactivate


###############################################################################
# STEP 5: Build a local pangenome graph with PGGB
#
# Parameters:
# -p 99     minimum identity for alignments
# -s 10000  segment length
# -c 10     graph sparsification
# -k 19     minimizer size
#
# Generates:
# Graph (.gfa/.og)
# VCF
# Various PGGB intermediate files
###############################################################################

conda activate pggb

pggb \
    -i tef_TT2_chr4B_merged.fa.gz \
    -o tef_TT2_chr4B \
    -n 26 \
    -t 2 \
    -p 99 \
    -s 10000 \
    -c 10 \
    -k 19 \
    -V 'Quncho_chr:1000'

###############################################################################
# STEP 6: Convert graph to ODGI format and simplify path names
###############################################################################

cd /home/shared/shared_all/e.teff/impg/tef_TT2_chr4B

conda deactivate
conda activate odgi

# List graph paths (accessions)
odgi paths \
    -i tef_TT2_chr4B_merged.fa.gz.f48ebc9.37c1a4f.ddda278.smooth.final.og \
    -L

# Export graph to GFA
odgi view \
    -i tef_TT2_chr4B_merged.fa.gz.f48ebc9.37c1a4f.ddda278.smooth.final.og \
    -g \
    > tef_TT2_chr4B_merged.gfa

# Simplify path names by removing suffixes after "_"
awk -F'\t' 'BEGIN{OFS=FS}{sub(/_.*/, "", $2); print}' \
    tef_TT2_chr4B_merged.gfa \
    > tef_TT2_chr4B_renamed.gfa

# Rebuild ODGI graph using renamed paths
odgi build \
    -g tef_TT2_chr4B_renamed.gfa \
    -o tef_TT2_chr4B_renamed.og

# Verify resulting path names
odgi paths -i tef_TT2_chr4B_renamed.og -L


###############################################################################
# STEP 7: Sort graph and generate graph visualization
###############################################################################

odgi sort \
    -i tef_TT2_chr4B_renamed.og \
    -o tef_TT2_chr4B_sorted.og \
    -O

odgi viz \
    -i tef_TT2_chr4B_sorted.og \
    -p accessions_order.txt \
    -F colors.tsv \
    -o tef_TT2_chr4B_sorted.viz.v1.png


###############################################################################
# STEP 8: Compute graph statistics
###############################################################################

# Obtain graph size and structural statistics
odgi stats \
    -i tef_TT2_chr4B_sorted.og \
    -S


###############################################################################
# STEP 9: Project gene and TE annotations onto graph coordinates
###############################################################################

odgi position \
    -i tef_TT2_chr4B_sorted.og \
    -b ../anno_impg.bed \
    -t 2

# Example mappings:
#
# Boni,7426,+    Boni,15623,+      Quncho_4B_003433
# Boni,9786,+    Boni,14923,+      Copia_LTR_retrotransposon


###############################################################################
# STEP 10: Extract subgraphs corresponding to annotation features
###############################################################################

# Gene subgraph
odgi extract \
    -i tef_TT2_chr4B_sorted.og \
    -r Boni:7426-15623 \
    -o Quncho_4B_003433.og

# TE subgraph
odgi extract \
    -i tef_TT2_chr4B_sorted.og \
    -r Boni:9786-14923 \
    -o Copia_LTR.og


###############################################################################
# STEP 11: Export feature subgraphs to GFA for Bandage
###############################################################################

odgi view \
    -i Quncho_4B_003433.og \
    -g \
    > Quncho_4B_003433.gfa

odgi view \
    -i Copia_LTR.og \
    -g \
    > Copia_LTR.gfa


###############################################################################
# STEP 12: Create Bandage annotation files
#
# CSV format:
# Node,Annotation
###############################################################################

grep '^S' Quncho_4B_003433.gfa |
cut -f2 |
awk 'BEGIN{print "Node,Annotation"} {print $1",Quncho_4B_003433"}' \
    > Quncho_4B_003433.csv

grep '^S' Copia_LTR.gfa |
cut -f2 |
awk 'BEGIN{print "Node,Annotation"} {print $1",Copia_LTR_retrotransposon"}' \
    > Copia_LTR.csv


###############################################################################
# STEP 13: Export complete graph to GFA
###############################################################################

odgi view \
    -i tef_TT2_chr4B_sorted.og \
    -g \
    > tef_TT2_chr4B_sorted.gfa


###############################################################################
# STEP 14: Generate a 2D graph layout
#
# anno_impg.bed is used to highlight genomic features on the graph
###############################################################################

odgi layout \
    -i tef_TT2_chr4B_sorted.og \
    -o tef_TT2_chr4B_sorted.lay

odgi draw \
    -i tef_TT2_chr4B_sorted.og \
    -c tef_TT2_chr4B_sorted.lay \
    -p tef_TT2_chr4B_sorted_2D.png \
    -C \
    -w 50 \
    -b ../anno_impg.bed