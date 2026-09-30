#!/usr/bin/env bash

#########################
#de novo genome assembly#
#########################
## assembly using hifiasm

#$geno.fq.gz --> raw reads
hifiasm -o $geno -t 28 --hg-size 0.6g -z 100 --primary $geno.fq.gz


# assembly quality evaluation
# BUSCO
busco -i geno.fa -c 10 -o geno -m genome -l poales_odb10

# QUAST
quast.py geno.fa -o geno_output

# LAI (After TE annotation)
LAI -genome geno.fasta.mod -intact geno.fasta.mod.EDTA.raw/LTR/geno.fasta.mod.pass.list -all geno.fasta.mod.EDTA.anno/geno.fasta.mod.out -t 20 -o LAI_geno


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


###################################
### Assembly_QC
#######################################

source ~/miniforge3/etc/profile.d/conda.sh
conda activate merqury

#echo "Starting Meryl: $(date)"

cd /mnt/nfs/home/riccucci/teff/assembly/merqury/out1
/mnt/nfs/home/riccucci/meryl-1.4.1/bin/meryl count k=21 threads=${SLURM_CPUS_PER_TASK} output quncho.illumina.meryl /mnt/nfs/home/riccucci/teff/assembly/quncho_shortreads/261668_ID3463_1-QUNCHO21_S178_L006_R1_001.fastq.gz /mnt/nfs/home/riccucci/teff/assembly/quncho_shortreads/261668_ID3463_1-QUNCHO21_S178_L006_R2_001.fastq.gz
#echo "Meryl completed: $(date)"

#echo "Starting Merqury: $(date)"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

/mnt/nfs/home/riccucci/merqury/merqury.sh quncho.illumina.meryl /mnt/nfs/home/riccucci/teff/assembly/Quncho_chr.filt.fasta quncho

#echo "Merqury completed: $(date)"
