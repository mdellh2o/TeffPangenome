README file

## Teff pangenome analysis scripts 
These scripts follow the analysis of the teff pangenome paper, using a mix of bash, R and python. Scripts follow are ordered in semi-chronological order following the flow of the paper. 
01_Genome assembly: de novo genome assembly, assembly quality evaluation, scaffolding
02_ Divergence Estimation: assembly Preparation and Per-Chromosome Divergence Estimation (Mash)
03_Genome annotation: gene annotation, repeat elements annotation, identification of centromeric and telomeric regions, filtering of genes overlapping TEs
04_Orthofinder: protein sequence extraction, orthogroup inference (OrthoFinder), core, dispensable and unique gene sets, functional annotation (eggNOG-mapper, InterProScan)
05_Functional annotation: functional annotation of core, dispensable and unique genes, using through InterProScan and EggNOG-mapper
06_GO_enrichment analysis: GO term enrichment of core, dispensable and unique genes
07_Ultrametric Species Tree: OrthoFinder and ultrametric species tree, preparation of the CAFE5 input files, ready-to-submit CAFE5 job script
08_Gene family analyses: gene family expansion and contraction analysis with CAFE5 (run and plots)
09_Gene expansion and contraction: significantly expanded and contracted gene lists from CAFE5, GO enrichment of expanded and contracted gene families
10_Structural variants: whole-genome alignment, structural variant calling (SyRI), plotting (plotsr)
11_RNA seq: RNA-seq alignment
12.DNA methylation analysis: ONT methylation BAM processing, DNA methylation analysis,  metagene DNA methylation and TE analysis for and  multi- genotype and  of core, dispensable and unique genes
13_Methylation gene expression: TE presence vs methylation, core vs dispensable gene expression, methylation-expression relationship
14_Haplotype analysis: k-mer identification based on IBSpy, by chromosome interval,  with phenotype and without phenotype
15_Combinatorial search: exhaustive combinatorial search to identify small haplotypic combinations with regards to traits of seed color and panicle shape
16_Classifiers:development of predictive models for target traits: seed color and panicle shape (27 pangenome accessions to 150 EIAR accessions) 
17_Machine learning method comparison: comparison of predictive models (nested cross-validation) for haplotype-based prediction of the phenotype
18_ Local pangenome graph:  A local pangenome graph for a single locus of the TT2 Locus (Chr 4B)  across 26 teff pangenome lines. 

