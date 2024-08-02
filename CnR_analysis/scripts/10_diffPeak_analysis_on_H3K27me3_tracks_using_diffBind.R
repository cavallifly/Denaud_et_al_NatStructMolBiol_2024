library(DiffBind)
library(profileplyr)
library(GenomicRanges)
library(ChIPseeker)
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(clusterProfiler)
library(ggplot2)
library("GenomicFeatures")
library("BRGenomics")
library(ReactomePA)
library(AnnotationDbi)
library(GenomicTools.fileHandler)
library(enrichplot)
library(biomaRt)
library(wordcloud)
library(edgeR)
library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(hrbrthemes)
library(ggpmisc)

comparisons <- c("Discs_WT_vs_double","Embryos_WT_vs_double")

for(comparison in comparisons)
{
    inDir <- paste0(comparison)
    setwd(inDir)

    samples <- read.delim("sample_table.tsv")
    K27me3 <- dba(sampleSheet = samples)

    #Counting reads in peaks, summits option is in case you want to center the counting in a range from the summit
    K27me3.counted <- dba.count(K27me3, summits=FALSE, filter=0) #no filter for peaks with low overlapping read counts, summit is for count in a given range from the summit
    K27me3.counted #FRiP: proportion of reads for that sample that overlap a peak in the consensus peakset

    K27me3.counted <- dba.contrast(K27me3.counted, contrast=c("Condition","Double","WT"))
    K27me3.counted
    dba.show(K27me3.counted,bContrasts=TRUE)
    K27me3.analyzed <- dba.analyze(K27me3.counted, method=DBA_DESEQ2, bBlacklist= FALSE, bGreylist = FALSE) #  DESEq2 without blacklist
    dba.show(K27me3.analyzed,bContrasts=TRUE) # see the amount of differential peaks

    #write report
    report_K27me3 <- dba.report(K27me3.analyzed, contrast = 1, method=DBA_DESEQ2, th=1, bCounts=TRUE)
    report_K27me3 <- as.data.frame(report_K27me3)  
    write.table(report_K27me3, "report_Embryo__H3K27me3_Double_vs_WT.tsv", sep="\t", quote=F, row.names=F)

    #scatter plot normalized counts
    plot_K27me3_double <- cbind(report_K27me3, dac_dom= 1) # add column for dac ploting color
    plot_K27me3_double$Domain<- with(plot_K27me3_double, paste0(seqnames, start, end)) #merge columns to id the domains
    plot_K27me3_double[plot_K27me3_double$Domain =="chr2L1635735016487149", "dac_dom"] <- 2


    colors <- c("gray", "blue")
    colors <- colors[as.numeric(plot_K27me3_double$dac_dom)]

    # Doing the scatter plot
    p3 <- ggplot(plot_K27me3_double, aes(x=Conc_WT, y=Conc_Double)) +
       geom_point(size=2) +
       geom_point(color=colors) +
       geom_abline(intercept = 0, slope = 1) +
       ylab("log2 Conc. H3K27me3 - Double") +
       xlab("log2 Conc. H3K27me3 - WT")
    pdf(file = paste0("scatterPlot_Conc_K27me3_",comparison,".pdf")   
    print(p3)
    dev.off()

    ## Doing the MAplot
    p3_MA <- ggplot(plot_K27me3_double, aes(x=Conc, y= Fold)) +
        geom_point(size=2,alpha=0.5,na.rm=T) +
	geom_point(color=colors) +
	geom_abline(intercept = 0, slope = 0) +
	geom_abline(intercept = -0.58, slope = 0, linetype = "dashed", color= "gray") +
	geom_abline(intercept = 0.58, slope = 0, linetype = "dashed", color= "gray") +
	ylim(-3,3) +
	ggtitle("H3K27me3 - Double vs WT") +
	ylab("log2 Fold change") +
	xlab("log2 Conc")
    pdf(file = "MAplot_Conc_K27me3_",comparison,".pdf")
    print(p3_MA)
    dev.off()
}