#!/usr/bin/env bash

#########################
#de novo genome assembly#
#########################
#hifiasm
#$geno.fq.gz --> raw reads
hifiasm -o $geno -t 28 --hg-size 0.6g -z 100 --primary $geno.fq.gz


#############
#Scaffolding#
#############
#ragtag

ragtag.py scaffold \
  -u \
  -t 32 \
  -q 40 \
  -f 25000 \
  -o /out_path \
# remove the unassigned contings and adjust headers 
grep -A 1 '>Quncho_' ../ragtag.scaffold.fasta > $geno.scaffolds.fasta
sed -i 's/_RagTag//' $geno.scaffolds.fasta
