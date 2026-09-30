## STAR Index generation

conda activate /projects/assembly_long_reads/conda_envs/parabricks_star_index/
mkdir -p /projects/assembly_long_reads/eragrostis_tef/quncho/quncho_star_index/

STAR \
--runMode genomeGenerate \
--genomeFastaFiles /projects/assembly_long_reads/eragrostis_tef/tef_varieties/assemblies/assemblies.filt/Salk_teff_dabbi_3.0.filt.fasta \
--runThreadN 20 \
--genomeDir /projects/assembly_long_reads/eragrostis_tef/dabbi/dabbi_star_index/ \
--genomeSAindexNbases 13



## Alignement step with Parabricks
# https://docs.nvidia.com/clara/parabricks/latest/ 

module load apptainer/1.3.4-gcc-14.2.0-rjsoecu
SIFFILES=/projects/novabreed/share/gmagris/software/sif_images/


# RNA igatech sequenced 
for id in /projects/assembly_long_reads/eragrostis_tef/quncho/quncho_RNA/RNA-Seq_data/*_R1_001.fastq.gz
do
    sample=$(basename $id _R1_001.fastq.gz)
    echo $sample 

    # align RNAseq data 
    singularity exec \
    -B /projects:/projects \
    --pwd /projects/assembly_long_reads/eragrostis_tef/quncho/ \
    --nv ${SIFFILES}/parabricks_4.5.0-1.sif pbrun rna_fq2bam \
    --in-fq /projects/assembly_long_reads/eragrostis_tef/quncho/quncho_RNA/RNA-Seq_data/${sample}_R1_001.fastq.gz /projects/assembly_long_reads/eragrostis_tef/quncho/quncho_RNA/RNA-Seq_data/${sample}_R2_001.fastq.gz \
    --genome-lib-dir /projects/assembly_long_reads/eragrostis_tef/quncho/quncho_star_index/ \
    --output-dir /projects/assembly_long_reads/eragrostis_tef/quncho/alignments/rna/ \
    --ref /projects/assembly_long_reads/eragrostis_tef/tef_varieties/assemblies/assemblies.filt/Quncho_chr.filt.fasta \
    --out-bam /projects/assembly_long_reads/eragrostis_tef/quncho/alignments/rna/${sample}.parabricks.bam \
    --read-files-command zcat    
done



## merge and filter fq.gz filles
conda activate /projects/assembly_long_reads/conda_envs/misc_tools/

samtools merge -@ 32 --write-index -o quncho_merged.parabricks.bam /projects/assembly_long_reads/eragrostis_tef/quncho/alignments/rna/177*.bam 
samtools view -@ 32 -e "endpos-pos<50000" -o quncho_merged_filt.parabricks.bam quncho_merged.parabricks.bam 
samtools index -@ 32 quncho_merged_filt.parabricks.bam

