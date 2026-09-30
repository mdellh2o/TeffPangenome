These scripts follow the analyses of the teff pangenome paper, using a mix of bash, R and Python. 
Scripts are ordered in semi-chronological order following the flow of the paper. 
Scripts are the following:

01_Genome assembly
de novo genome assembly

02_Genome annotation
Gene annotation
Repeart Elements annotation
Identification of centromeric and telomeric regions

03_Orthofinder 
Protein sequence extraction
Functional annotation: eggNOG-mapper
Functional annotation: InterProScan

04_ Functional Annotation
Functional annotation of core, dispensable and unique genes

05_GO_enrichment analysis
GO term enrichment of core, dispensable and unique genes

06_Ultrametric Species Tree
OrthoFinder and ultrametric species tree
Preparation of the CAFE5 input files
Ready-to-submit CAFE5 job script
CAFE run

07_Gene Family Analysis
Gene family expansion and contraction analysis with CAFE5

08_Gene expansion and contraction
Significantly expanded and contracted gene lists from CAFE5

09_Structural Variants
Whole-genome alignment
Structural variant calling SyRI

10_RNA_seq
RNA-seq Alignment 

11_DNA methylation analysis 
Extract gene sets from OrthoFinder
DNA methylation meta-analysis
Transposable element meta-analysis
Identification of core, dispensable and unique genes
Metagene DNA methylation analysis
Metagene transposable element analysis

12_Methylation_Gene_expression
TE presence vs methylation. 
Core vs dispensable gene expression. 
Methylation-expression relationship.

13_Haplotype analysis
kmer identification based on IBSpy, by chromosome interval

14_Combinatorial search
Observation of haplotypic combinations with regards to traits

15_Classifiers
Developement of predictive model(s) for target traits

16_Machine learning method comparison
Comparison of predictive model(s)
