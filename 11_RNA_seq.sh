#!/usr/bin/env bash
###################################
#RNA-seq Alignment
###############################

## STAR Index generation

conda activate /PATH/parabricks_star_index/
mkdir -p /PATH/quncho/quncho_star_index/

STAR \
--runMode genomeGenerate \
--genomeFastaFiles /PATH/Salk_teff_dabbi_3.0.filt.fasta \
--runThreadN 20 \
--genomeDir /PATH/dabbi_star_index/ \
--genomeSAindexNbases 13



## Alignement step with Parabricks
# https://docs.nvidia.com/clara/parabricks/latest/ 

module load apptainer/1.3.4-gcc-14.2.0-rjsoecu
SIFFILES=/PATH/sif_images/


# RNA igatech sequenced 
for id in /PATH/RNA-Seq_data/*_R1_001.fastq.gz
do
    sample=$(basename $id _R1_001.fastq.gz)
    echo $sample 

    # align RNAseq data 
    singularity exec \
    -B /projects:/projects \
    --pwd /PATH/quncho/ \
    --nv ${SIFFILES}/parabricks_4.5.0-1.sif pbrun rna_fq2bam \
    --in-fq /PATH/quncho_RNA/RNA-Seq_data/${sample}_R1_001.fastq.gz /PATH/RNA-Seq_data/${sample}_R2_001.fastq.gz \
    --genome-lib-dir /PATH/quncho_star_index/ \
    --output-dir /PATH/alignments/rna/ \
    --ref /PATH/Quncho_chr.filt.fasta \
    --out-bam /PATH/rna/${sample}.parabricks.bam \
    --read-files-command zcat    
done

## merge and filter fq.gz filles
conda activate /projects/assembly_long_reads/conda_envs/misc_tools/

samtools merge -@ 32 --write-index -o quncho_merged.parabricks.bam /PATH/177*.bam 
samtools view -@ 32 -e "endpos-pos<50000" -o quncho_merged_filt.parabricks.bam quncho_merged.parabricks.bam 
samtools index -@ 32 quncho_merged_filt.parabricks.bam

