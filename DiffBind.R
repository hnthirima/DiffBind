source("https://bioconductor.org/biocLite.R")
biocLite("DiffBind")

library(DiffBind)

#1 Setting the path
path = setwd('/home/hthirima/hthirima_fast/2018_CR/030918_DiffBind_H3K27ac/DiffBind3/')

#2 Reading sample sheet to a dataframe
samples <- read.csv("SampleSheet_DiffBind_030918_H3K27ac_MACS2_nrrwpeak.csv")
names(samples)
samples

#3 Peaksets are read directly using the samplesheet
H3K27ac_030918 <- dba(sampleSheet="SampleSheet_DiffBind_030918_H3K27ac_MACS2_nrrwpeak.csv")

#4 Tells you how many peaks are identified in each peakset
#In the first line: total number of unique peaks after merging overlapping ones (in paranthesis)
#Dimenstions of binding matrix (#of samples vs sites overlapped in at least 2 samples)
H3K27ac_030918

#5 Correlation heatmap using occupancy (peak caller score) data
plot(H3K27ac_030918)

#6 Counting reads
H3K27ac_030918.count <- dba.count(H3K27ac_030918)

#7 The number of consensus peaks is shown in the first line. 
#FRiP (Fraction of Reads In Peaks): indicate which samples show more enrichemnt overall
H3K27ac_030918.count

#8 Correlation heatmap using affinity (read count) data
plot(H3K27ac_030918.count)

#Don't do step 9 and 10 if you want all pairwise constrasts to be considered.

#9 Establishing contrast. We tell DiffBind which samples fall in which groups
H3K27ac_030918_contrast <- dba.contrast(H3K27ac_030918.count, categories = DBA_CONDITION)

#10 Perfroming differential analysis. 
#The main differential analysis function is invoked as follows. 
#This will run an DESeq2 analysis using the default matrix (FDR <= 0.05)
H3K27ac_030918_analyze <- dba.analyze(H3K27ac_030918_contrast)
H3K27ac_030918_analyze

# Alternative to step 9 and 10
H3K27ac_030918_analyze <- dba.analyze(H3K27ac_030918.count)

#11 Correlation heatmap based on analysis
plot(H3K27ac_030918_analyze)

#12 Retreiving differentially bound sites
H3K27ac_030918.DB <- dba.report(H3K27ac_030918_analyze)

#Extract diff bound sites for each contrast
WT.R.DB <- dba.report(H3K27ac_030918_analyze, contrast = 1)
WT.NL.DB <- dba.report(H3K27ac_030918_analyze, contrast = 2)
R.NL.DB <- dba.report(H3K27ac_030918_analyze, contrast = 3)

WT.NL.DB
WT_NL_th1.DB = dba.report(WT_NL.DB, th=9.9E-05)

WT_NL.DB <- dba.report(WT_NL.DB)
WT_NL.DB.Gain <- WT_NL.DB(WT_NL.DB$Fold>0, th=9.9E-05)

WT_NL_th1.DB.Gain <- dba.report(H3K27ac_030918_analyze, contrast = 2, th = 9.9E-05, bGain = TRUE)

write.csv(H3K27ac_030918.DB, file ="DiffBoundSites_H3K27ac_030918_2.csv")
write.csv(WT.NL.DB, file = "DiffBoundSites_H3K27ac_WTNL_030918_DB3.csv")

write.csv(WT_NL_th1.DB.Gain, file = "DB2_H3K27ac_WTNL_th1_Gain.csv")

#13 Plotting
#MA plot: use data(XXXXX_analysis)
dba.plotMA(H3K27ac_030918_analyze, contrast = 3)
dba.plotMA(H3K27ac_030918_analyze, bXY=TRUE)

#14 PCA plot
dba.plotPCA(H3K27ac_030918.count,DBA_TISSUE,label=DBA_CONDITION)

dba.plotPCA(H3K27ac_030918_analyze, label=DBA_TISSUE)

dba.plotPCA(H3K27ac_030918_analyze, attributes=c(DBA_TISSUE,DBA_CONDITION), label= DBA_REPLICATE)

#15 Box plot
sum(H3K27ac_030918.DB$Fold<0)
sum(H3K27ac_030918.DB$Fold>0)

pvals <- dba.plotBox(H3K27ac_030918.DB)
pvals

#16 Heatmaps
corvals <- dba.plotHeatmap(H3K27ac_030918_analyze, correlations = FALSE)

