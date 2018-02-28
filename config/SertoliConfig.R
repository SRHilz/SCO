#Config for Sertoli-enriched analysis
rawCountsFile <- 'data/merged_counts_noNC.txt'
specificListFile <- 'data/20180209_SertoliEnriched_HumanGeneSymbols.txt'
tag <- '20180223_SertoliSubset'
ofInterest <- list(Transmembrane=c("CLDN11","PVRL2","GJB1","CDH2","JAM3","TUBG2"), Adapter=c("TJP3","PKP1","PKP3","AXIN1"))
toFilter <- 1