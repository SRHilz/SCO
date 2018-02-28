```{r include=FALSE, cache=FALSE}
library(edgeR)
library(goseq)
library(ggplot2)
library(dplyr)
```

```{r include=FALSE, cache=FALSE}
## FUNCTIONS
sem <- function(x) sd(x)/sqrt(length(x))

plotPCA <- function(CPMOutputFile, outfile){
  par(bg = 'white')
  par(mfrow=c(1,1))
  par(mar=c(4,4,4,4))
  data=read.table(CPMOutputFile,header = TRUE, sep = ',')
  names = gsub('TS0','',colnames(data)[2:length(colnames(data))])
  filtered_data <- data[which(rowSums(data[2:length(colnames(data))]>1) >= 4),]
  nrm_count_matrix=as.matrix(filtered_data[,2:length(colnames(data))])
  log=log10(1+nrm_count_matrix)
  tlog=t(log)
  pca=prcomp(tlog)
  summary(pca)
  colors = c( "dodgerblue2", 
              "dodgerblue2",  
              "dodgerblue2", 
              "dodgerblue2", 
              "firebrick", 
              "firebrick",
              "firebrick",
              "firebrick",
              "firebrick", 
              "firebrick",
              "firebrick")
  raw <- pca$x[,1:3]
  xaxis <- 1
  yaxis <- 2
  xlab <- paste('PC',toString(xaxis),
                ' (',toString(round(100*summary(pca)$importance[2,xaxis])),'%)')
  ylab <- paste('PC',toString(yaxis),
                ' (',toString(round(100*summary(pca)$importance[2,yaxis])),'%)')
  xaxisExtra <-(max(raw[,xaxis])-min(raw[,xaxis]))/4
  yaxisExtra <-(max(raw[,yaxis])-min(raw[,yaxis]))/4
  print(plot(raw[,xaxis], raw[,yaxis],  col=colors, pch=20, 
       xlim=c(min(raw[,xaxis])-xaxisExtra,max(raw[,xaxis])+xaxisExtra),
       ylim=c(min(raw[,yaxis])-yaxisExtra,max(raw[,yaxis])+yaxisExtra),
       xlab=xlab,
       ylab=ylab,
       cex=2,
       cex.lab=1,
       cex.axis=1))
  print(text(raw[,xaxis], raw[,yaxis],names, c(1,1), pos=2))
  print(legend(-3,2,c("Normal/Hypo","SCO"),pch=c(20,20),cex=1.3,
               col=c("dodgerblue2","firebrick"),bty='n'))
  dev.copy2pdf(file=outfile)
}

plotBarchart <- function(geneList, group, cpms, outfile){
  geneList <- geneList %>% unlist
  subset <- cpms[geneList,]
  meansNormal <- rowMeans(subset[,which(group==1)])
  meansSCO <- rowMeans(subset[,which(group==2)])
  semNormal <- apply(subset[,which(group==1)],1,sem)
  semSCO <- apply(subset[,which(group==2)],1,sem)
  toPlotNormal <- cbind(meansNormal,semNormal,names(meansNormal),
                        rep('Normal',length(meansNormal)))
  toPlotNormal <- as.data.frame(toPlotNormal)
  colnames(toPlotNormal) <- c('MeanCPM','SEM','Gene','Group')
  toPlotSCO <- cbind(meansSCO,semSCO,names(meansSCO),rep('SCO',length(meansSCO)))
  toPlotSCO <- as.data.frame(toPlotSCO)
  colnames(toPlotSCO) <- c('MeanCPM','SEM','Gene','Group')
  toPlot <- rbind(toPlotNormal, toPlotSCO)
  toPlot$MeanCPM <- as.numeric(as.character(toPlot$MeanCPM))
  toPlot$SEM <- as.numeric(as.character(toPlot$SEM))
  rownames(toPlot) <- NULL
  toPlot$Gene <- factor(toPlot$Gene, levels = geneList)
  print(ggplot(data=toPlot, aes(x=Gene, y=MeanCPM, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=MeanCPM-SEM, ymax=MeanCPM+SEM), width=.1, 
                  position=position_dodge(.9)) +
    scale_y_continuous(trans="log2") +
    scale_fill_manual(values=c("dodgerblue2", "firebrick")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20, color='black'), 
          axis.title = element_text(size = 20), 
          axis.text.y = element_text(size=20, color='black'), 
          legend.text=element_text(size=20), 
          panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.title=element_blank()))
  dev.copy2pdf(file=outfile)
}

plotScatter <- function(cpms, group, outfile){
  meansNormal <- log(rowMeans(cpms[,which(group==1)]),2)
  meansSCO <- log(rowMeans(cpms[,which(group==2)]),2)
  plot(x=meansNormal, y=meansSCO, pch=1, cex=.5, 
       xlab="log2 CPM Normal",ylab="log2 CPM SCO")
  abline(lm(meansSCO~meansNormal))
  rsq <- round(summary(lm(meansSCO~meansNormal))$r.squared,4)
  legend("topleft",bty='n',legend=paste0("Rsquared = ",rsq))
  dev.copy2pdf(file=outfile)
}

annotateGO <- function(GO.wall,invertedGetgo){
  GO.wall$foreground_genes <- 'NA'
  for (cat in GO.wall$category){
    GO.wall[which(GO.wall$category==cat),]$foreground_genes <- paste(invertedGetgo[[cat]],collapse=';')
  }
  return(GO.wall)
}
  
invertGetgo <- function(getgolist){
  inversion <- list()
  getgolist <- getgolist[!is.na(names(getgolist))]
  for (gene in names(getgolist)){
    for (goterm in getgolist[[gene]]){
      if (!goterm %in% names(inversion)){
        inversion[[goterm]] = c()
      }
      inversion[[goterm]] <- append(inversion[[goterm]], gene)
    }
  }
  return(inversion)
}

performGO <- function(binaryList, outfile){
  print("Table of input values")
  print(table(binaryList))
  pwf=nullp(binaryList,'hg19',"geneSymbol")
  GO.wall=goseq(pwf,"hg19","geneSymbol")
  print("Top 20 most significant GO terms")
  top <- GO.wall[,c(6,2)]
  colnames(top) <- c("term","pvalue")
  rownames(top) <- NULL
  print(head(top,20))
  getgolist <- getgo(names(binaryList[which(binaryList==1)]), 'hg19','geneSymbol')
  getGeneList <- invertGetgo(getgolist)
  GO.wall.anno <- annotateGO(GO.wall,getGeneList)
  write.csv(GO.wall.anno,outfile, quote=F)
}
```

## MAIN
```{r message=FALSE}
# read in config file for analysis - change config to analyze a different
#  subset of genes (choices for config - "Sertoli", "Leydig", "Union")
config <- "Union"
source(paste0('config/',
              config,'Config.R'))

# set up file names
CPMOutputFile <- paste0('output/',tag,'_CPMs.csv')
DAOutputFile <- paste0('output/',tag,'_NHvSCO_edgeR_results.csv')
PCAOutputFile <- paste0('output/',tag, '_PCA.pdf')
MAFile <- paste0('output/',tag, '_logCPM_v_logFC.pdf')
BoxplotRawOutputFile <- paste0('output/',tag, '_PreNorm_CPMs.pdf')
BoxplotNormOutputFile <- paste0('output/',tag, '_PostNorm_CPMs.pdf')
GOSpecificFile <- paste0('output/',tag, '_SpecificSubset_GeneOntology.csv')
GOUpFile <- paste0('output/',tag, '_Up_GeneOntology.csv')
GODownFile <- paste0('output/',tag, '_Down_GeneOntology.csv')
scatterFile <- paste0('output/',tag, '_Scatter.pdf')

# read in raw counts file
sampleTable_edgeR<-read.delim(rawCountsFile, row.names='gene')

# check dimensions
dim(sampleTable_edgeR)

# build logical vector of rownames that are not genes but summary outputs of HTSeq
noint = rownames(sampleTable_edgeR) %in% c("__ambiguous",
                                           "__too_low_aQual", 
                                           "__not_aligned", 
                                           "__no_feature",
                                           "__alignment_not_unique")

# set grouping - first four are normal, remaining are SCO
group<-factor(c(1,1,1,1,2,2,2,2,2,2,2))

# build DGEList object
d<-DGEList(counts=sampleTable_edgeR,group=group)

# subset original matrix by genes that are expressed over a CPM cutoff, and, 
  # if toFilter==1, that are in the provided gene list
if (toFilter==1){
  specific_list <- scan(file=specificListFile, what=character())
  specific = toupper(rownames(sampleTable_edgeR)) %in% toupper(specific_list)
  paste0('In specific list: ',
        length(specific_list[toupper(specific_list) %in% toupper(rownames(sampleTable_edgeR))]))
  paste('Not in specific list: ',
        length(specific_list[!toupper(specific_list) %in% toupper(rownames(sampleTable_edgeR))]))
  keep <- !noint & specific 
  }else{
  keep <- !noint
}
d<- d[keep,]

# check dimensions after filtering
dim(d)

# perform GO on specific gene list compared to all genes
if (toFilter == 1){
  specificGenes=as.integer(rownames(sampleTable_edgeR) %in% specific_list)
  names(specificGenes) <- rownames(sampleTable_edgeR)
  performGO(specificGenes, GOSpecificFile)
}

# look at CPMs of selected genes before any normalization
boxplot(log(cpm(d)+.001,2), col='red', ylab="log2 CPM")
dev.copy2pdf(file=BoxplotRawOutputFile)

# calculate the normalization factors (this will correct for overall differences in count 
  # means between samples)
d<-calcNormFactors(d, method="RLE")#normalizing by log median

# show the normalization factor calculated for each library
d$samples

# look at CPMs of selected genes after normalization to mean counts
boxplot(log(cpm(d)+.001,2), col='red', ylab="log2 CPM")
dev.copy2pdf(file=BoxplotNormOutputFile)

# estimate common dispersion across all samples
d<-estimateCommonDisp(d)

# view common dispersion
sqrt(d$common.disp)

# estimate individual dispersion for each gene
d<-estimateTagwiseDisp(d)

# ouptut CPMs to file
write.csv(cpm(d),CPMOutputFile)

# examine CPMs as normal vs SCO scatterplot
plotScatter(cpm(d), group, scatterFile)

# examine CPMs for proteins of interest
for (list in names(ofInterest)){
  print(paste0("Generating barplot for ",list," proteins"))
  barFile <- paste0("output/",tag,"_",list,"_Bar.pdf")
  plotBarchart(ofInterest[list], group, cpm(d), barFile)
}

# PubQuality PCA plot with % explained
plotPCA(CPMOutputFile, PCAOutputFile)

# default of exactTest uses tag dispresion, does pairwise comp, 
  # comparing 2 to 1 Normal + Hypo vs SCO
NHvSCO_edgeR=exactTest(d, pair=c("1","2"))

# format results
results_NHvSCO<-topTags(NHvSCO_edgeR, n = nrow( NHvSCO_edgeR$table ) )$table

# make a vector of all differentially expressed genes
NHvSCO_detags <- rownames(results_NHvSCO)[results_NHvSCO$FDR < 0.05]

# summarize results
summary(decideTestsDGE(NHvSCO_edgeR, p=0.05, adjust="BH"))

# make a MA style plot
plotSmear(NHvSCO_edgeR, de.tags=NHvSCO_detags)
dev.copy2pdf(file=MAFile)

# output to a file
write.csv(results_NHvSCO,DAOutputFile)

# perform GO on significantly upregulated genes
genes_up_NHvSCO=as.integer(results_NHvSCO$logFC > 0 & 
                             rownames(results_NHvSCO) %in% NHvSCO_detags)
names(genes_up_NHvSCO) <- rownames(results_NHvSCO)
performGO(genes_up_NHvSCO, GOUpFile)

# perform GO on significantly downregulated genes
genes_down_NHvSCO=as.integer(results_NHvSCO$logFC < 0 & 
                               rownames(results_NHvSCO) %in% NHvSCO_detags)
names(genes_down_NHvSCO) <- rownames(results_NHvSCO)
performGO(genes_down_NHvSCO, GODownFile)
```

