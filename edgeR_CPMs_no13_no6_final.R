source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("goseq")
biocLite("geneLenDataBase")
biocLite("org.Hs.eg.db")
biocLite("mygene")
biocLite("biomaRt")
library(edgeR)#call edgeR
library(goseq)
library(geneLenDataBase)
library(org.Hs.eg.db)
library(mygene)
library(biomaRt)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("pheatmap")) {
  install.packages("pheatmap", dependencies = TRUE)
  library(pheatmap)
}
library(grid)

plotPCACorr <- function(x, y, xlab, ylab){
  corr <- cor(x, y, method = "pearson")
  print(corr)
  plot(x,y, xlab=xlab, ylab=ylab)
  print(c(par('usr')))
  position <- c(par("usr")[2]-(par("usr")[2]/3),par("usr")[4]-(par("usr")[4]/3))
  print(position)
  print(round(corr,2))
  text(position[1],position[2],paste('r=',toString(round(corr,2))))
}

#Config for Sertoli-specific analysis
rawCountsFile <- 'merged_counts_noNC_wCLU.txt'
specificListFile <- '20170101_SertoliSpecific_HumanGeneSymbols_Revised.txt'
CPMOutputFile <- '20170101_RNA_CPM_Values_noNCnew_No13_No6_wCLU_SertoliNormed_SertoliSubset_Revised.csv'
DAOutputFile <- "NHvSCO_edgeR_results_No13_No6_wCLU_SertoliNormed_SertoliSubset_20170101.csv"
toFilter <- 1

#Config for Leydig-specific analysis
rawCountsFile <- 'merged_counts_noNC_wCLU.txt'
specificListFile <- '20170101_LeydigSpecific_HumanGeneSymbols_Revised.txt'
CPMOutputFile <- '20170101_RNA_CPM_Values_noNCnew_No13_No6_wCLU_LeydigNormed_LeydigSubset_Revised.csv'
DAOutputFile <- "NHvSCO_edgeR_results_No13_No6_wCLU_LeydigNormed_LeydigSubset_20170101.csv"
toFilter <- 1

#Config for non-specific analysis
rawCountsFile <- 'merged_counts_noNC_wCLU.txt'
CPMOutputFile <- '20170101_RNA_CPM_Values_noNCnew_No13_No6_wCLU.csv'
DAOutputFile <- "NHvSCO_edgeR_results_No13_No6_wCLU_20170101.csv"
toFilter <- 0

#Main
sampleTable_edgeR<-read.delim(rawCountsFile, row.names='gene')
dim(sampleTable_edgeR)
noint = rownames(sampleTable_edgeR) %in% c("__ambiguous","__too_low_aQual", "__not_aligned", "__no_feature","__alignment_not_unique")
sampleTable_edgeR = sampleTable_edgeR[1:12]#use this if you want to remove the 13th sample, the tech rep (can keep if want to assess technical reporducibility)
sampleTable_edgeR$TS006=NULL#and use this to remove TS006, the sample that was unable to be placed in a binary categorization of "SCO" and "nonSCO"
group<-factor(c(1,1,1,1,2,2,2,2,2,2,2))
d<-DGEList(counts=sampleTable_edgeR,group=group)
if (toFilter==1){
  specific_list <- scan(file=specificListFile, what=character())
  specific = toupper(rownames(sampleTable_edgeR)) %in% toupper(specific_list)
  paste('In specific list: ',length(specific_list[toupper(specific_list) %in% toupper(rownames(sampleTable_edgeR))]),sep='')
  paste('Not in specific list: ',length(specific_list[!toupper(specific_list) %in% toupper(rownames(sampleTable_edgeR))]),sep='')
  keep <- rowSums(cpm(d)>.4) >= 2 & !noint & specific #filters out anything without CPM>n in at least two libraries, our minimal rep size. .4 is a CPM cutoff for raw counts of ~6 for a ~20 million read library; also removes pesky rows that are not genes (like __ambiguous)
}else{
  keep <- rowSums(cpm(d)>.4) >= 2 & !noint #filters out anything without CPM>n in at least two libraries, our minimal rep size. .4 is a CPM cutoff for raw counts of ~6 for a ~20 million read library; also removes pesky rows that are not genes (like __ambiguous)
}
d<- d[keep,]
dim(d)
d<-calcNormFactors(d)#normalizing by dfault
d<-estimateCommonDisp(d)
d<-estimateTagwiseDisp(d)
sqrt(d$common.disp)
plotMDS(d)
write.csv(cpm(d),CPMOutputFile)

#PubQuality PCA plot with % explained
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
colors = c( "seagreen", "seagreen",  "seagreen", "seagreen", "firebrick", "firebrick","firebrick","firebrick","firebrick", "firebrick","firebrick")
raw <- pca$x[,1:3]
xaxis <- 1
yaxis <- 2
xlab <- paste('PC',toString(xaxis),' (',toString(round(100*summary(pca)$importance[2,xaxis])),'%)')
ylab <- paste('PC',toString(yaxis),' (',toString(round(100*summary(pca)$importance[2,yaxis])),'%)')
xaxisExtra <-(max(raw[,xaxis])-min(raw[,xaxis]))/4
yaxisExtra <-(max(raw[,yaxis])-min(raw[,yaxis]))/4
plot(raw[,xaxis], raw[,yaxis],  col=colors, pch=20, 
       xlim=c(min(raw[,xaxis])-xaxisExtra,max(raw[,xaxis])+xaxisExtra),
       ylim=c(min(raw[,yaxis])-yaxisExtra,max(raw[,yaxis])+yaxisExtra),
       xlab=xlab,
       ylab=ylab,
       cex=2,
       cex.lab=1,
       cex.axis=1)
text(raw[,xaxis], raw[,yaxis],names, c(1,1))
#legend(-3,2,c("Normal/Hypo","SCO"),pch=c(20,20),cex=1.3,col=c("seagreen","firebrick"),bty='n')

#Boxplot CPMs
data=read.table(CPMOutputFile,header = TRUE, sep = ',')
samples = colnames(data)[2:length(colnames(data))]
row.names(data)=data$X
data=data[,2:length(colnames(data))]
boxplot(data)# not too informative, because tons of highly expressed outliers
boxplot(data, ylim=c(0,70), col='red')#zoom in on boxplots, ignoring highly expressed outliers
boxplot(data+.001, ylim=c(0,70), col='red')#here I add a psuedocount, just so I can then take a log below; this shows this psuedocount obviously has no effect on the data trend
boxplot(log(data+.001),2, col='red')

#DA analysis - Normal + Hypo vs SCO
group<-factor(c(1,1,1,1,2,2,2,2,2,2,2))
NHvSCO_edgeR=exactTest(d, pair=c("1","2"))#default of exactTest uses tag dispresion, does pairwise comp, comparing 2 to 1
NHvSCO_detags <- rownames(topTags(NHvSCO_edgeR, n=20))
cpm(d)[NHvSCO_detags,]#couples with the above line, shows 
summary(NHvSCO_de <- decideTestsDGE(NHvSCO_edgeR, p=0.05, adjust="BH"))
NHvSCO_detags <- rownames(d)[as.logical(NHvSCO_de)]
plotSmear(NHvSCO_edgeR, de.tags=NHvSCO_detags)
results_NHvSCO<-topTags(NHvSCO_edgeR, n = nrow( NHvSCO_edgeR$table ) )$table
write.csv(results_NHvSCO,DAOutputFile)
genes_NHvSCO=as.integer(p.adjust(NHvSCO_edgeR$table$PValue[NHvSCO_edgeR$table$logFC!=0],method="BH")<.05)

#Make Scatter of  specific normalized
data=read.table('20170101_RNA_CPM_Values_noNCnew_No13_No6_wCLU_SertoliNormed_SertoliSubset_Revised.csv',header = TRUE, sep = ',', row.names="X")
sig = rownames(data) %in% NHvSCO_detags
normal_data = data[,1:4]
sco_data = data[,5:11]
normal_sig_data = rowSums(normal_data[sig,])/ncol(normal_data[sig,])
normal_notsig_data = rowSums(normal_data[!sig,])/ncol(normal_data[!sig,])
sco_sig_data = rowSums(sco_data[sig,])/ncol(sco_data[sig,])
sco_notsig_data = rowSums(sco_data[!sig,])/ncol(sco_data[!sig,])
plot(log(normal_notsig_data,2),log(sco_notsig_data,2),ylim=c(-7,12),xlim=c(-7,12), xlab="Normal/Hypo log2 CPM",ylab="SCO log2 CPM")
points(log(normal_sig_data,2),log(sco_sig_data,2),col="red")
#legend(0,10,c("Not Significantly Different","Significantly Differend"),pch=c(1,1),cex=1,col=c("black","red"),bty='n')

#Make Scatter of non-sertoli specific normalized
data=read.table('20160524_RNA_CPM_Values_noNCnew_No13_No6_wCLU.csv',header = TRUE, sep = ',', row.names="X")
sig = rownames(data) %in% NHvSCO_detags
normal_data = data[,1:4]
sco_data = data[,5:11]
normal_sig_data = rowSums(normal_data[sig,])/ncol(normal_data[sig,])
normal_notsig_data = rowSums(normal_data[!sig,])/ncol(normal_data[!sig,])
sco_sig_data = rowSums(sco_data[sig,])/ncol(sco_data[sig,])
sco_notsig_data = rowSums(sco_data[!sig,])/ncol(sco_data[!sig,])
plot(log(normal_notsig_data,2),log(sco_notsig_data,2),ylim=c(-7,12),xlim=c(-7,12), xlab="Normal/Hypo log2 CPM",ylab="SCO log2 CPM")
points(log(normal_sig_data,2),log(sco_sig_data,2),col="red")
#legend(0,10,c("Not Significantly Different","Significantly Differend"),pch=c(1,1),cex=1,col=c("black","red"),bty='n')
