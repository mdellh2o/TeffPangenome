
#!/usr/bin/env bash

##Protein sequence extraction## 

#Protein sequences were extracted for each genotype from the EDTA-filtered Helixer gene annotation. 
#Chromosome FASTA files are first rewrapped to 60 characters per line with seqkit (required by AGAT), 
#then AGAT extracts the protein sequences into ${geno}.protein.fa.

#agat needs 60 characters per line

for file in /PATH/fasta/*.fasta
do
base=$(basename "$file")
geno=${base%%_*}
#geno=$(basename $file _chr.filt.fasta)
fa_dir=/PATH/fasta/
prot_dir=/PATH/protein_sequences/
anno_dir=/PATH/helixer_annotation_EDTAfilt/

conda activate seqkit
cd ${fa_dir}

seqkit seq -w 60 ${base} > ${geno}.60w.fasta
conda deactivate
conda activate agat

cd ${prot_dir}
agat_sp_extract_sequences.pl \
-g ${anno_dir}${geno}.filt.EDTA.gff3 \
-f ${fa_dir}${geno}.60w.fasta \
-p \
-o ${geno}.protein.fa
done



#OrthoFinder

##Orthogroups were inferred from the protein sequences of all genotypes with OrthoFinder (DIAMOND search, MAFFT alignments and IQ-TREE gene trees). 

#!/usr/bin/env bash

conda deactivate
conda activate orthofinder

orthofinder -f ${prot_dir} -M msa -S diamond -t 96 -a 8 -I 1.5 -T iqtree -A mafft


orthofinder \
 -f /PATH/orthofinder \
 -X \
 -t 29


##Core, dispensable and unique gene sets##

#All protein files are merged, and the proteins of the core, dispensable and unique genes 
#are extracted with seqkit using the gene ID lists. 


cd /PATH/protein_fasta_files
cat *.fa > all_proteins.fa


conda activate seqkit

seqkit grep -f /PATH/core_genes.txt all_protein.fa > core_genes.fa
seqkit grep -f /PATH/disp_genes.txt all_protein.fa > disp_genes.fa
seqkit grep -f /PATH/uniq_genes.txt all_protein.fa > uniq_genes.fa


cd /PATH/protein_fasta_files/

sed -i 's/*//' core_genes.fa
sed -i 's/*//' disp_genes.fa
sed -i 's/*//' uniq_genes.fa

##Functional annotation: InterProScan##

#Protein domains and GO terms were annotated for each gene set with InterProScan 

for set in core_genes disp_genes uniq_genes
do
cd /PATH/functional_anno
#cd functional_anno/v0/interpro/my_interproscan/interproscan-5.77-108.0/interproscan.sh
#InterProScan
sh v0/interpro/my_interproscan/interproscan-5.77-108.0/interproscan.sh -i /PATH/protein_sequences/protein_fasta_files/${set}.fa -f tsv -goterms --cpu 28 -o v1/interpro.${set}.txt # --> output: interpro.core_genes.txt, interpro.dispensable_genes.txt, interpro.unique_genes.txt
done

##Functional annotation: eggNOG-mapper##

#Each gene set was also annotated with eggNOG-mapper

#########server embeds
source ~/miniforge3/etc/profile.d/conda.sh
conda create -n eggnog

for set in core_genes disp_genes uniq_genes
do
cd /PATH/functional_anno
#Eggnog
emapper.py --data_dir ./eggnog-mapper-data -i ./${set}.fa -o ${set}.eggnog --cpu 28 # --> output: core_genes.eggnog.emapper.annotations, dispensable_genes.eggnog.emapper.annotations, unique_genes.eggnog.emapper.annotations
done

##protein sequences for each orthogroup##
#retrieve protein sequences for each OG
#tentativo di recuperare tutte le sequenze... troppo lungo, Giorgia dovrebbe passarti direttamente gli output di orthofinder

#for file in /PATH/OG/*.txt
#do
#conda activate seqkit
#og=$(basename "$file" .txt)
#prot_dir=/PATH/protein_sequences/protein_fasta_files
#og_dir=/PATH/OG

#echo ${prot_dir}/all_protein.fa
#echo ${og_dir}/${og}.fa
#done
#seqkit grep -f ${file} ${prot_dir}/all_protein.fa > ${og_dir}/fasta/${og}.fa
#done

#############
#quindi hai deciso di selezionare randomicamente un singolo gene per ciascun ortogruppo
#--> vedi 'check_core_disp_uniqu.R'
for file in /PATH/*in_og.txt
do
conda activate seqkit
og=$(basename "$file" .txt)
prot_dir=/PATH/protein_sequences/protein_fasta_files
og_dir=/PATH/OG

seqkit grep -f ${file} ${prot_dir}/all_protein.fa > ${og_dir}/fasta/${og}.fa
done


##Functional annotation of OG representatives#

#InterProScan and eggNOG-mapper are run only on the representative gene of each orthogroup 
#(core, dispensable; unique is commented out).

#invece di fare interproscan e eggnog per tutti i geni fallo per i rappresentanti di ogni OG

sed -i 's/*//' single_core_in_og.fa
sed -i 's/*//' single_disp_in_og.fa
sed -i 's/*//' single_uniq_in_og.fa

for set in single_core_in_og single_disp_in_og #single_uniq_in_og
do
conda deactivate
fasta_dir=/PATH/OG/fasta
intepro_dir=/PATH/interpro/my_interproscan/interproscan-5.77-108.0
out_dir_interpro=/PATH/interpro_out
eggnog_dir=/PATH/eggnog-mapper-2.1.14
out_dir_eggnog=/PATH/eggnog_out

#InterProScan
sh ${intepro_dir}/interproscan.sh -i ${fasta_dir}/${set}.fa -f tsv -goterms --cpu 28 -o ${out_dir_interpro}/interpro.${set}.txt # --> output: interpro.core_genes.txt, interpro.dispensable_genes.txt, interpro.unique_genes.txt

#Eggnog
conda activate eggnog
${eggnog_dir}/emapper.py --data_dir ${eggnog_dir}/data -i ${fasta_dir}/${set}.fa -o ${out_dir_eggnog}/${set}.eggnog --cpu 28 # --> output: core_genes.eggnog.emapper.annotations, dispensable_genes.eggnog.emapper.annotations, unique_genes.eggnog.emapper.annotations

done