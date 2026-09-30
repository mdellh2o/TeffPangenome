###-------------------------------###
# NEW CORE DISPENSABLE UNIQUE GENES #
###-------------------------------###

# based on orthogroups_all.tsv from orthofinder, extract the genes in each category (core, dispensable, unique) for Quncho

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



###-----------------------###
# METHYLATION META-ANALYSIS #
###-----------------------###

cd /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure

for gff in /home/shared/shared_all/e.teff/genes_anno/helixer_annotation_EDTAfilt/*.filt.EDTA.gff3
do
base=$(basename "$gff")
genos="${base%%.*}"
# Convert gff3 in bed file for positive and negative strand separately, subgenosma A and B separately
awk -F'\t' -v OFS='\t' -v genos="$genos" '

BEGIN {
    corefile = "/home/shared/shared_all/e.teff/functional_anno/v1/ms_figure/" genos "_core.txt"
    dispfile = "/home/shared/shared_all/e.teff/functional_anno/v1/ms_figure/" genos "_dispensable.txt"
    uniqfile = "/home/shared/shared_all/e.teff/functional_anno/v1/ms_figure/" genos "_unique.txt"

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


# methylation bed file
cd /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure

for METH in /home/shared/shared_all/e.teff/methylation_self_alignments/*.pileup.bed
#METH=/home/shared/shared_all/e.teff/methylation_self_alignments/*.pileup.bed
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

##################
#sei arrivato qua# --> ricordati di mettere i ${genos} negli output
##################


# TE bed file
#from gff to bed

cd /home/shared/shared_all/e.teff/TE_anno/EDTA_all_genosme/EDTA_IGA

for gff in *_chr.filt.fasta.mod.EDTA.TEanno.gff3
do
base=$(basename "$gff")
genos="${base%%.*}"

awk 'BEGIN{OFS="\t"}
     !/^#/ {
         print $1, $4-1, $5, $9, ".", $7
     }' ${gff} > ${genos}_TE.bed
done

cd /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure
for TE in /home/shared/shared_all/e.teff/TE_anno/EDTA_all_genosme/EDTA_IGA/*_TE.bed
do
#TE=/home/shared/shared_all/e.teff/TE_anno/EDTA_all_genosme/EDTA_IGA
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
cd /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure
for TE in /home/shared/shared_all/e.teff/TE_anno/EDTA_all_genosme/EDTA_IGA/*_TE.bed
do
#TE=/home/shared/shared_all/e.teff/TE_anno/EDTA_all_genosme/EDTA_IGA
for tbed in /home/shared/shared_all/e.teff/functional_anno/v1/ms_figure/*core_A_positive_transcript_with_flanks_bins.bed
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

## Plotting in R

#conda activate /projects/assembly_long_reads/conda_envs/epigenosmics/

R
setwd("/home/shared/shared_all/e.teff/functional_anno/v1/ms_figure/")
library(ggplot2)
library(dplyr)
library(viridis)

#methylation data

read_set <- function(prefix, subgen){

  bind_rows(lapply(c("core","dispensable","unique"), function(cat){

    f <- paste0(cat,"_",subgen,"_CG_percentile_meth.txt")

    df <- read.table(f, stringsAsFactors = FALSE)
    colnames(df) <- c("region","bin","meth")

    df$category <- cat
    df$subgenosme <- subgen

    df
  }))
}

df <- bind_rows(
  read_set(NULL,"A"),
  read_set(NULL,"B")
)

df$x <- with(df,
  ifelse(region == "a",
         bin - 100,        # -100 → 0 (TSS)
  ifelse(region == "b",
         bin,              # 0 → 100 (gene body)
  ifelse(region == "c",
         bin + 100, NA)))  # 100 → 200 (TTS downstream)
)

df$category <- factor(df$category,
                      levels=c("core","dispensable","unique"))
df$subgenosme <- factor(df$subgenosme,
                       levels = c("A","B"),
                       labels = c("Subgenosme A","Subgenosme B"))
df$meth <- df$meth / 100

p <- ggplot(df, aes(x = x, y = meth, color = category)) +

  geom_line(linewidth = 1, alpha = 0.9) +

  facet_wrap(~subgenosme, nrow = 1) +

  scale_color_viridis_d(option = "magma", begin = 0.1, end = 0.95) +

  scale_x_continuous(
    breaks = c(-100, 0, 50, 100, 200),
    labels = c("-2kb", "TSS", "gene body",
               "TTS", "+2kb")
  ) +

  theme_classic() +

  labs(
    x = "",
    y = "DNA methylation level (mC/C)",
    color = "Gene category"
  ) +

  theme(
    strip.text = element_text(face = "bold", size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

ggsave("Quncho_meta_gene_A_B.png",
       p, width = 12, height = 6, dpi = 300)
	   
####################################################################################

library(ggplot2)
library(dplyr)
library(viridis)

# Read one category/subgenosme combination
read_set <- function(geno, subgen){

  bind_rows(
    lapply(c("core","dispensable","unique"), function(cat){

      f <- paste0(geno, "_", cat, "_", subgen,
                  "_CG_percentile_meth.txt")

      df <- read.table(f, stringsAsFactors = FALSE)
      colnames(df) <- c("region","bin","meth")

      df$category <- cat
      df$subgenosme <- subgen

      df
    })
  )
}

# Find genos automatically
files <- list.files(pattern = "_CG_percentile_meth\\.txt$")

genos <- unique(
  sub("_(core|dispensable|unique)_[AB]_CG_percentile_meth\\.txt$",
      "",
      files)
)

# Skip T33_correct
genos <- setdiff(genos, "T33_correct")
# Optional: remove incomplete file left in folder
genos <- genos[genos != "core_A_CG_percentile_meth.txt"]

print(genos)

# Generate one figure per geno
for(geno in genos){

  cat("Processing:", geno, "\n")

  df <- bind_rows(
    read_set(geno, "A"),
    read_set(geno, "B")
  )

  df$x <- with(df,
               ifelse(region == "a",
                      bin - 100,
                      ifelse(region == "b",
                             bin,
                             ifelse(region == "c",
                                    bin + 100, NA))))

  df$category <- factor(
    df$category,
    levels = c("core","dispensable","unique")
  )

  df$subgenosme <- factor(
    df$subgenosme,
    levels = c("A","B"),
    labels = c("Subgenosme A","Subgenosme B")
  )

  df$meth <- df$meth / 100

  p <- ggplot(df,
              aes(x = x,
                  y = meth,
                  colour = category,
                  group = category)) +
    geom_line(linewidth = 1, alpha = 0.9) +
    facet_wrap(~subgenosme, nrow = 1) +
    scale_color_viridis_d(
      option = "magma",
      begin = 0.1,
      end = 0.95
    ) +
    scale_x_continuous(
      breaks = c(-100,0,50,100,200),
      labels = c("-2kb","TSS",
                 "gene body","TTS","+2kb")
    ) +
    theme_classic() +
    labs(
      title = geno,
      x = "",
      y = "DNA methylation level (mC/C)",
      color = "Gene category"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 15),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 15),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    )

  ggsave(
    paste0(geno, "_meta_gene_A_B.png"),
    p,
    width = 12,
    height = 6,
    dpi = 300
  )
}

#######################################################################################
# TE data 

library(ggplot2)
library(dplyr)
library(viridis)

read_set <- function(geno, subgen){

  bind_rows(
    lapply(c("core", "dispensable", "unique"), function(cat){

      f <- paste0(
        geno, "_", cat, "_", subgen,
        "_TE_percentile.txt"
      )

      df <- read.table(f, stringsAsFactors = FALSE)
      colnames(df) <- c("region", "bin", "TE")

      df$category <- cat
      df$subgenome <- subgen

      df
    })
  )
}

# Detect genos automatically
files <- list.files(pattern = "_TE_percentile\\.txt$")

genos <- unique(
  sub(
    "_(core|dispensable|unique)_[AB]_TE_percentile\\.txt$",
    "",
    files
  )
)

# Skip T33_correct
genos <- setdiff(genos, "T33_correct")

for(geno in genos){

  cat("Processing:", geno, "\n")

  df <- bind_rows(
    read_set(geno, "A"),
    read_set(geno, "B")
  )

  df$x <- with(
    df,
    ifelse(region == "a",
           bin - 100,
           ifelse(region == "b",
                  bin,
                  ifelse(region == "c",
                         bin + 100, NA)))
  )

  df$category <- factor(
    df$category,
    levels = c("core", "dispensable", "unique")
  )

  df$subgenome <- factor(
    df$subgenome,
    levels = c("A", "B"),
    labels = c("Subgenome A", "Subgenome B")
  )

  p <- ggplot(
    df,
    aes(x = x, y = TE, color = category)
  ) +
    geom_line(linewidth = 1, alpha = 0.9) +
    facet_wrap(~subgenome, nrow = 1) +
    scale_color_viridis_d(
      option = "magma",
      begin = 0.1,
      end = 0.95
    ) +
    scale_x_continuous(
      breaks = c(-100, 0, 50, 100, 200),
      labels = c(
        "-2kb", "TSS",
        "gene body", "TTS", "+2kb"
      )
    ) +
    theme_classic() +
    labs(
      title = geno,
      x = "",
      y = "TE coverage (%)",
      color = "Gene category"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 15),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 15),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    )

  ggsave(
    paste0(geno, "_meta_gene_TE_A_B.png"),
    p,
    width = 12,
    height = 6,
    dpi = 300
  )
}

####################################################################################################


##make unique figure across genomes
library(ggplot2)
library(dplyr)
library(viridis)

read_set <- function(geno, subgen){

  bind_rows(
    lapply(c("core","dispensable"), function(cat){

      f <- paste0(
        geno, "_", cat, "_", subgen,
        "_TE_percentile.txt"
      )

      df <- read.table(f, stringsAsFactors = FALSE)
      colnames(df) <- c("region","bin","TE")

      df$category <- cat
      df$subgenome <- subgen
      df$geno <- geno

      df
    })
  )
}

files <- list.files(pattern = "_TE_percentile\\.txt$")

genos <- unique(
  sub(
    "_(core|dispensable|unique)_[AB]_TE_percentile\\.txt$",
    "",
    files
  )
)

# remove T33_correct
genos <- setdiff(genos, "T33_correct")

# read all genos
df <- bind_rows(
  lapply(genos, function(geno){

    bind_rows(
      read_set(geno, "A"),
      read_set(geno, "B")
    )

  })
)

# metagene coordinate
df$x <- with(
  df,
  ifelse(region == "a",
         bin - 100,
         ifelse(region == "b",
                bin,
                ifelse(region == "c",
                       bin + 100, NA)))
)

# average across genos
df_mean <- df %>%
  group_by(subgenome, category, x) %>%
  summarise(
    TE = mean(TE, na.rm = TRUE),
    .groups = "drop"
  )

df_mean$category <- factor(
  df_mean$category,
  levels = c("core","dispensable")
)

df_mean$subgenome <- factor(
  df_mean$subgenome,
  levels = c("A","B"),
  labels = c("Subgenome A","Subgenome B")
)

p <- ggplot(
  df_mean,
  aes(x = x, y = TE, colour = category)
) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~subgenome, nrow = 1) +
  scale_color_viridis_d(
    option = "magma",
    begin = 0.1,
    end = 0.95
  ) +
  scale_x_continuous(
    breaks = c(-100,0,50,100,200),
    labels = c("-2kb","TSS","gene body","TTS","+2kb")
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Mean TE coverage (%)",
    color = "Gene category"
  )

ggsave(
  "Mean_TE_profile_all_genos.png",
  p,
  width = 12,
  height = 6,
  dpi = 300
)

saveRDS(df_mean, "Mean_TEcoverage_all_genos_df_mean.rds")
saveRDS(p, "Mean_TEcoverage_all_genos_plot.rds")

#####METHYLATION

read_set <- function(geno, subgen){

  bind_rows(
    lapply(c("core","dispensable"), function(cat){

      f <- paste0(
        geno, "_", cat, "_", subgen,
        "_CG_percentile_meth.txt"
      )

      df <- read.table(f, stringsAsFactors = FALSE)
      colnames(df) <- c("region","bin","meth")

      df$category <- cat
      df$subgenome <- subgen
      df$geno <- geno

      df
    })
  )
}

files <- list.files(pattern = "_CG_percentile_meth\\.txt$")

genos <- unique(
  sub(
    "_(core|dispensable|unique)_[AB]_CG_percentile_meth\\.txt$",
    "",
    files
  )
)

# remove T33_correct
genos <- setdiff(genos, "T33_correct")

# read all genos
df <- bind_rows(
  lapply(genos, function(geno){

    bind_rows(
      read_set(geno, "A"),
      read_set(geno, "B")
    )

  })
)

# metagene coordinate
df$x <- with(
  df,
  ifelse(region == "a",
         bin - 100,
         ifelse(region == "b",
                bin,
                ifelse(region == "c",
                       bin + 100, NA)))
)

# average across genos
df_mean <- df %>%
  group_by(subgenome, category, x) %>%
  summarise(
    meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

df_mean$meth <- df_mean$meth / 100

df_mean$subgenome <- factor(
  df_mean$subgenome,
  levels = c("A","B"),
  labels = c("Subgenome A","Subgenome B")
)

df_mean$category <- factor(
  df_mean$category,
  levels = c("core", "dispensable")
)

p <- ggplot(
  df_mean,
  aes(x = x, y = meth, colour = category)
) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~subgenome, nrow = 1) +
  scale_color_viridis_d(
    option = "magma",
    begin = 0.1,
    end = 0.95
  ) +
  scale_x_continuous(
    breaks = c(-100,0,50,100,200),
    labels = c("-2kb","TSS","gene body","TTS","+2kb")
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "DNA methylation level (mC/C)",
    color = "Gene category"
  )

ggsave(
  "Mean_DNA_meth_all_genos.png",
  p,
  width = 12,
  height = 6,
  dpi = 300
)

saveRDS(df_mean, "Mean_DNA_meth_all_genos_df_mean.rds")
saveRDS(p, "Mean_DNA_meth_all_genos_plot.rds")