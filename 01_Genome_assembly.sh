
#!/usr/bin/env bash

#########################
#de novo genome assembly#
#########################
## assembly using hifiasm

#geno.fq.gz --> raw reads
hifiasm -o $geno -t 28 --hg-size 0.6g -z 100 --primary $geno.fq.gz

##Assembly quality evaluation##

# assembly quality evaluation
# BUSCO
busco -i geno.fa -c 10 -o geno -m genome -l poales_odb10

# QUAST
quast.py geno.fa -o geno_output

# LAI (After TE annotation)
LAI -genome geno.fasta.mod -intact geno.fasta.mod.EDTA.raw/LTR/geno.fasta.mod.pass.list -all geno.fasta.mod.EDTA.anno/geno.fasta.mod.out -t 20 -o LAI_geno

##Scaffolding

#Contigs were ordered and oriented into chromosome-level scaffolds with RagTag.

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

##Assembly QC (Merqury)##

#Assembly consensus quality (QV) and k-mer completeness were evaluated with Merqury. 
#A k-mer database (k=21) is first built from the Illumina short reads with Meryl

###################################
### Assembly_QC
#######################################

source ~/miniforge3/etc/profile.d/conda.sh
conda activate merqury

#echo "Starting Meryl: $(date)"

cd PATH/merqury/out1
/PATH/meryl-1.4.1/bin/meryl count k=21 threads=${SLURM_CPUS_PER_TASK} output quncho.illumina.meryl /PATH/261668_ID3463_1-QUNCHO21_S178_L006_R1_001.fastq.gz /PATH/quncho_shortreads/261668_ID3463_1-QUNCHO21_S178_L006_R2_001.fastq.gz
#echo "Meryl completed: $(date)"

#echo "Starting Merqury: $(date)"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

/PATH//merqury/merqury.sh quncho.illumina.meryl /PATH/Quncho_chr.filt.fasta quncho

#echo "Merqury completed: $(date)"


### Estimation of assemblies divergence


#!/usr/bin/env bash
/PATH/chr_genome
#header should follow this scheme:
#[sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]

ls *_chr.filt.fasta | while read f; do
    sample_name=$(echo $f | cut -f 1 -d '.');
    echo ${sample_name}
    fastix -p "${sample_name}#1#" $f >> E_tef_pangenome.fasta
done
bgzip -@ 4 E_tef_pangenome.fasta
samtools faidx E_tef_pangenome.fasta.gz

#Assuming we will work with one chromosome at a time, we estimate the sequence divergence for each set of chromosomes. To partition the sequences by chromosome, execute:
cut -f 1 E_tef_pangenome.fasta.gz.fai | cut -f 3 -d '#' | sort | uniq | while read CHROM; do
    CHR_FASTA=E_tef_pangenome.$CHROM.fasta.gz
    samtools faidx E_tef_pangenome.fasta.gz $(grep -P $"$CHROM\t" E_tef_pangenome.fasta.gz.fai | cut -f 1) | bgzip -@ 4 > $CHR_FASTA

    echo "Generated $CHR_FASTA"
done

#to estimate the distance of each input sequence to every other sequence in the set, we use mash. In particular, the mash triangle command outputs a lower-triangular distance matrix. 
for chr in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B 8A 8B 9A 9B 10A 10B
do
mash triangle E_tef_pangenome.$chr.fasta.gz > E_tef_pangenome.$chr.mash_triangle.txt
cat E_tef_pangenome.$chr.mash_triangle.txt | column -t
done


for chr in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B 8A 8B 9A 9B 10A 10B
do
sed -i 's/#1#//' E_tef_pangenome.$chr.mash_triangle.txt
done





 

