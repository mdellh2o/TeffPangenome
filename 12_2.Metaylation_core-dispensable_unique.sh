#!/usr/bin/env bash

## Core, dispensable and unique genes per genotype

###-------------------------------###
# NEW CORE DISPENSABLE UNIQUE GENES #
###-------------------------------###

#based on orthogroups_all.tsv from orthofinder, extract the genes in each category (core, dispensable, unique) for Quncho
#Metagene profiles of CG DNA methylation and TE coverage
#(2 kb upstream, gene body, 2 kb downstream) for core, dispensable and unique genes, in subgenome A and B of each genotype.

cd /home/shared/shared_all/e.teff/functional_anno/v1
#quncho col 6 --> $6


mkdir -p ms_figure

for col in {2..28}
do
    colname=$(awk -F'\t' -v col="$col" '
        NR==1 {
            x=$col
            gsub(/[^[:alnum:]_.-]/, "_", x)
            print x
            exit
        }
    ' orthogroups_all.tsv)

    awk -F'\t' \
        -v col="$col" \
        -v colname="$colname" '
    BEGIN {
        while (getline x < "core_orthogroups.txt")        core[x]=1
        while (getline x < "dispensable_orthogroups.txt") disp[x]=1
        while (getline x < "unique_orthogroups.txt")      unique[x]=1
    }

    NR > 1 {
        if ($(col) == "") next

        n = split($(col), g, ", ")

        for (i=1; i<=n; i++) {
            if ($1 in core)
                print g[i] >> "ms_figure/" colname "_core.txt"

            if ($1 in disp)
                print g[i] >> "ms_figure/" colname "_dispensable.txt"

            if ($1 in unique)
                print g[i] >> "ms_figure/" colname "_unique.txt"
        }
    }
    ' orthogroups_all.tsv
done

## Gene BED files with flanking regions
## Converts the mRNA features of each genotype to BED files, separately for category, subgenome (A or B) and strand. Each gene gets three regions: 
## upstream 2 kb (a), gene body (b) and downstream 2 kb (c), oriented by strand


###-----------------------###
# METHYLATION META-ANALYSIS #
###-----------------------###

cd /PATH/ms_figure

for gff in /PATH/*.filt.EDTA.gff3
do
base=$(basename "$gff")
genos="${base%%.*}"
# Convert gff3 in bed file for positive and negative strand separately, subgenosma A and B separately
awk -F'\t' -v OFS='\t' -v genos="$genos" '

BEGIN {
    corefile = "/PATH/ms_figure/" genos "_core.txt"
    dispfile = "/PATH/ms_figure/" genos "_dispensable.txt"
    uniqfile = "/PATH/ms_figure/" genos "_unique.txt"

    print "corefile:", corefile > "/dev/stderr"

    while ((getline line < corefile) > 0) {
        class[line] = "core"
        ncore++
    }

    while ((getline line < dispfile) > 0) {
        class[line] = "dispensable"
        ndisp++
    }

    while ((getline line < uniqfile) > 0) {
        class[line] = "unique"
        nuniq++
    }

    print genos, ncore, ndisp, nuniq > "/dev/stderr"
}

$3=="mRNA" {

    if(match($9,/ID=([^;]+)/,m))
        gene=m[1]
    else next

    sub(/\..*$/,"",gene)

    if(!(gene in class))
        next

    if(gene ~ /_[0-9]+A_/)
        subgen="A"
    else if(gene ~ /_[0-9]+B_/)
        subgen="B"
    else
        next

    strand=($7=="+" ? "positive" : "negative")

    outfile=genos"_"class[gene]"_"subgen"_"strand"_transcript_with_flanks.bed"
	print outfile > "/dev/stderr"
    chr=$1
    start=$4-1
    end=$5

    print chr,start,end,"b",$7 >> outfile

    if($7=="+"){
        a_start=(start-2000<0 ? 0 : start-2000)
        a_end=start-1

        c_start=end+1
        c_end=end+2000
    }
    else{
        a_start=end+1
        a_end=end+2000

        c_start=(start-2000<0 ? 0 : start-2000)
        c_end=start-1
    }

    print chr,a_start,a_end,"a",$7 >> outfile
    print chr,c_start,c_end,"c",$7 >> outfile
}
' "${gff}"
done


## Metagene bins
## Divides each region into 100 bins with bedtools
# make windows of 100 bins for each strand separately
conda activate bedtools
for bed in *_transcript_with_flanks.bed
do
    if [[ $bed == *negative* ]]; then
        bedtools makewindows -n 100 -i srcwinnum -reverse -b "$bed" \
        | sed 's/_/\t/1' > "${bed%.bed}_bins.bed"
    else
        bedtools makewindows -n 100 -i srcwinnum -b "$bed" \
        | sed 's/_/\t/1' > "${bed%.bed}_bins.bed"
    fi
done

## DNA methylation per bin
# methylation bed file
cd /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure

for METH in /PATH/*.pileup.bed
#METH=/PATH/*.pileup.bed
do
base=$(basename "$METH")
genos="${base%%.*}"

for class in core dispensable unique
do
    for subgen in A B
    do

        cat \
            ${genos}_${class}_${subgen}_positive_transcript_with_flanks_bins.bed \
            ${genos}_${class}_${subgen}_negative_transcript_with_flanks_bins.bed \
        | sort -k4,4 -k5,5n \
        > tmp_bins.bed

        bedtools intersect \
            -a tmp_bins.bed \
            -b ${METH} \
            -wa -wb \
        | sort -k4,4 -k5,5n \
        | bedtools groupby \
            -g 4,5 \
            -c 16 \
            -o mean \
        > ${genos}_${class}_${subgen}_CG_percentile_meth.txt

    done
done
done

rm tmp_bins.bed

## TE coverage per bin

##################
#sei arrivato qua# --> ricordati di mettere i ${genos} negli output
##################


# TE bed file
#from gff to bed

cd /PATH/EDTA_IGA

for gff in *_chr.filt.fasta.mod.EDTA.TEanno.gff3
do
base=$(basename "$gff")
genos="${base%%.*}"

awk 'BEGIN{OFS="\t"}
     !/^#/ {
         print $1, $4-1, $5, $9, ".", $7
     }' ${gff} > ${genos}_TE.bed
done

cd /PATH/ms_figure
for TE in /PATH/*_TE.bed
do
#TE=/PATH/EDTA_IGA
base=$(basename "$TE")
genos="${base%%_*}"
for class in core dispensable unique
do
    for subgen in A B
    do

        cat \
            ${genos}_${class}_${subgen}_positive_transcript_with_flanks_bins.bed \
            ${genos}_${class}_${subgen}_negative_transcript_with_flanks_bins.bed \
        | sort -k4,4 -k5,5n \
        > tmp_bins.bed

        bedtools intersect \
         -a tmp_bins.bed \
         -b "$TE" \
         -wao \
        | awk 'BEGIN{OFS="\t"} {len=$3-$2; ov=$NF; print $4,$5,len,ov}' \
        | sort -k1,1 -k2,2n \
        | bedtools groupby -g 1,2 -c 3,4 -o sum,sum \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,($4/$3)*100}' \
        > ${genos}_${class}_${subgen}_TE_percentile.txt

    done
done
done

###########v1
cd /PATH/ms_figure
for TE in /PATH/*_TE.bed
do
#TE=/PATH/EDTA_IGA
for tbed in /PATH/ms_figure/*core_A_positive_transcript_with_flanks_bins.bed
do
base=$(basename "$tbed")
genos="${base%%_*}"
for class in core dispensable unique
do
    for subgen in A B
    do

        cat \
            ${genos}_${class}_${subgen}_positive_transcript_with_flanks_bins.bed \
            ${genos}_${class}_${subgen}_negative_transcript_with_flanks_bins.bed \
        | sort -k4,4 -k5,5n \
        > tmp_bins_v1.bed

        bedtools intersect \
         -a tmp_bins_v1.bed \
         -b "$TE" \
         -wao \
        | awk 'BEGIN{OFS="\t"} {len=$3-$2; ov=$NF; print $4,$5,len,ov}' \
        | sort -k1,1 -k2,2n \
        | bedtools groupby -g 1,2 -c 3,4 -o sum,sum \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,($4/$3)*100}' \
        > ${genos}_${class}_${subgen}_TE_percentile_v1.txt

    done
done
done
done

rm tmp_bins.bed
