#Config for Sertoli-Leydig-specific Union analysis
rawCountsFile <- 'data/merged_counts_noNC.txt'
specificListFile <- 'data/20180209_UnionSertoliLeydig_HumanGeneSymbols.txt'
tag <- '20180223_UnionSubset'
ofInterest <- list(Transmembrane=c("CLDN11","PVRL2","GJB1","CDH2","JAM3","TUBG2"), Adapter=c("TJP3","PKP1","PKP3","AXIN1"), AndrogenBiosyn=c("STAR","CYP11A1","CYP17A1","HSD3B2","HSD17B3"))
toFilter <- 1
