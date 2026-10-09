
#!/usr/bin/env bash

### Estimation of assemblies divergence

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





 

