#!/usr/bin/env bash

#####################
#Structural variants#
#####################

minimap2 -a -x asm5 -t 22 --eqx Quncho.fa ${geno}.scaffolds.fasta > ${geno}.onquncho.sam
samtools view -@ 5 -Sb -o ${geno}.onquncho.bam ${geno}.onquncho.sam
samtools sort ${geno}.onquncho.bam -o ${geno}.onquncho_sorted.bam
samtools index ${geno}.onquncho_sorted.bam

############
#Syri
python /path/syri -c ${geno}.onquncho_sorted.bam -r Quncho.fasta -q ${geno}.scaffolds.fasta -F B --prefix ${geno}.
python /path/plotsr ${geno}.syri.out /Quncho.fa ${geno}.scaffolds.fasta -H 8 -W 5