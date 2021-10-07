source("https://bioconductor.org/biocLite.R")
biocLite("DiffBind")

library(DiffBind)

#1 Setting the path
path = setwd('/home/hthirima/hthirima_fast/091117_H3K4me/DiffBind/')

#2 Reading sample sheet to a dataframe
samples <- read.csv("H3K4me1and2_WTvsNL.csv")
names(samples)
samples

#3 Peaksets are read directly using the samplesheet
H3K4me1and2_WTvsNL <- dba(sampleSheet="H3K4me1and2_WTvsNL.csv")

#4 Tells you how intervals identified in each peakset
#In the first line: total number of unique peaks after merging overlapping ones
#Dimensions of binding matrix (#of samples vs sites overlap in at least 2 samples)
H3K4me1and2_WTvsNL

#5 Correlation heatmap using occupancy (peak caller score) data
plot(H3K4me1and2_WTvsNL)

#6 Counting reads
H3K4me1and2_WTvsNL <- dba.count(H3K4me1and2_WTvsNL)

#7 The number of consensus peaks set, FRiP
H3K4me1and2_WTvsNL

#8 Correlation heatmap using affinity (read count) data
plot(H3K4me1and2_WTvsNL)

#9 Establishing contrast. Input which samples fall in which groups
H3K4me1and2_WTvsNL <- dba.contrast(H3K4me1and2_WTvsNL, categories = DBA_CONDITION, minMembers = 2)

#10 Perfroming differential analysis. 
#The main differential analysis function is invoked as follows:
H3K4me1and2_WTvsNL <- dba.analyze(H3K4me1and2_WTvsNL)
H3K4me1and2_WTvsNL

#11 Correlation heatmap based on analysis
plot(H3K4me1and2_WTvsNL, contrast=1)

#12 Retreiving differentially bound sites
H3K4me1and2_WTvsNL.DB <- dba.report(H3K4me1and2_WTvsNL)

H3K4me1and2_WTvsNL.DB

write.csv(H3K4me1and2_WTvsNL.DB, file ="DiffBoundSites_me1.csv")

#13 Plotting
#MA plot: use data(XXXXX_analysis)
dba.plotMA(H3K4me1and2_WTvsNL, contrast = 1)

dba.plotMA(H3K4me1and2_WTvsNL, bXY=TRUE)

#14 PCA plot
dba.plotPCA(H3K4me1and2_WTvsNL, contrast = 1, label=DBA_FACTOR)

dba.plotPCA(H3K4me1and2_WTvsNL, attributes=c(DBA_CONDITION), label= DBA_FACTOR)

#15 Box plot
sum(H3K4me1and2_WTvsNL.DB$Fold<0)
sum(H3K4me1and2_WTvsNL.DB$Fold>0)

pvals <- dba.plotBox(H3K4me1and2_WTvsNL)
pvals

#16 Heatmaps
corvals <- dba.plotHeatmap(H3K4me1and2_WTvsNL, contrast=1, correlations = FALSE, scale="row")



