#!/usr/bin/env bash

#################
#Gene annotation#
#################

Helixer.py \
--fasta-path ${file} \
--lineage land_plant \
--gff-output-path ${out}.gff3 \
--temporary-dir ../gene_prediction/ \
--subsequence-length 64152 

################
#TEs annotation#
################

.conda/envs/EDTA/bin/EDTA.pl \
  --genome ${geno}.scaffolds.fasta \
  --species others \
  --step all \
  --cds path/Eragrostis_tef.Salk_teff_dabbi_3.0.cds.all.fa \
  --curatedlib path/12870_2016_725_MOESM1_ESM.fas \
  --anno 1 \
  --sensitive 1 \
  --evaluate 1 \
  --threads 28 \
  --overwrite 1

######################################################
#identification of centromeric and telomeric regions##
######################################################

#Tandem Repeat Finder
trf ${geno}.scaffolds.fasta 2 5 7 80 10 50 2000
#RepeatMasker
RepeatMasker --no_is -nolow -q -pa 26 --dir /path/ ${geno}.scaffolds.fasta
#Telomeric regions --> Tidk
tidk explore -w 100000 -c Poales ${geno}.scaffolds.fasta


##########################################################
##########################################################
## filter Helixer gene predictions based on TE annotation#
##########################################################
# removing TE-overlapping genes
# merge TE

for TE in ../*.gff3
do
    name=$(basename "$TE")
    name=${name%%[-_.]*}
    
    sort -k1,1 -k4,4n ${TE} > ${name}.EDTA.TEanno.sorted.gff3
    bedtools merge -i ${name}.EDTA.TEanno.sorted.gff3 > ${name}.EDTA.TEanno.merged.bed
done

for gene in ../*.filt.gff3
do
    name=$(basename "$gene")
    name=${name%%[-_.]*}
    TE=/projects/assembly_long_reads/eragrostis_tef/tef_varieties/TE_EDTA/${name}.EDTA.TEanno.merged.bed # TE 

    # identify genes that overlap with TEs more than 80% of their CDS length
    bedtools intersect -a ${gene} -b ${TE} -wao \
        > tmp/gene_TE_intersect.gff3

    awk -F'\t' '
    $3 == "CDS" {

        match($9, /ID=([^;]+)/, arr)
        cds_id = arr[1]

        gene_id = cds_id
        sub(/\.(CDS).*$/, "", gene_id)
        sub(/\.1$/, "", gene_id)

        # 1) sum all CDS overlap (also duplications)
        overlap_sum[gene_id] += $NF

        # 2) sum all CDS length (unique CDS only)
        cds_key = gene_id "|" cds_id
        if (!(cds_key in cds_seen)) {
            cds_length[gene_id] += ($5 - $4 + 1)
            cds_seen[cds_key] = 1
        }
    }
    END {
        for (g in overlap_sum) {
            ratio = (cds_length[g] > 0) ? overlap_sum[g] / cds_length[g] : 0
            if (ratio >= 0.8)
                print g
        }
    }
    ' tmp/gene_TE_intersect.gff3 \
    > tmp/ids_to_remove.txt

    # remove from the gff3
    awk -F'\t' '
    BEGIN {
        while ((getline line < "tmp/ids_to_remove.txt") > 0) {
            remove[line] = 1
        }
        close("tmp/ids_to_remove.txt")
    }
    {
        match($9, /ID=([^;]+)/, arr)
        id_full = arr[1]
        sub(/\..*$/, "", id_full)
        if (!(id_full in remove)) print $0
    }
    ' OFS="\t" ${gene} > ${name}.filt.EDTA.gff3

done
