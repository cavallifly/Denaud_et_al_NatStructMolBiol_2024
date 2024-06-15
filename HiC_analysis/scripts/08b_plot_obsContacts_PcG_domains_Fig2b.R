require(reshape2) # cat tabulate to matrix
library(ggplot2)
library(ggpubr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
source('./scripts/auxFunctions.R')
options("scipen"=999)
colors <- rainbow(10)

obsData <- list(
	LE_WT        = c("hic.hic_LE_WT_Rep1_dm6_YO.hic_LE_WT_Rep1_dm6_YO_1",
	   	         "hic.hic_LE_WT_Rep2_dm6_YO.hic_LE_WT_Rep2_dm6_YO_1",
			 "hic.hic_LE_WT_Rep3_dm6_VL.hic_LE_WT_Rep3_dm6_VL_1"),
	larvae_DWT   = c("hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1",
			 "hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1",
			 "hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1"),
        pupae_DWT    = c("hic.hic_pupae_DWT_Rep1_dm6_BeS.hic_pupae_DWT_Rep1_dm6_BeS_1",
			 "hic.hic_pupae_DWT_Rep2_dm6_BeS.hic_pupae_DWT_Rep2_dm6_BeS_1",
			 "hic.hic_pupae_DWT_Rep3_dm6_BeS.hic_pupae_DWT_Rep3_dm6_BeS_1")
        )
samples = names(obsData)
print(paste0("Samples ",samples))

refSample <- "LE_WT"
refChrom  <- "chr2L"

DomainsFile <- "./scripts/list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed"
Domains     <- read.table(DomainsFile, header=F)
colnames(Domains) <- c("chrom","start","end","domain")
Domains     <- Domains[Domains$chrom == refChrom,]
print(head(Domains))

RegionsFile <- "./scripts/list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed"
Regions     <- read.table(RegionsFile, header=F)
colnames(Regions) <- c("chrom","start","end","domain")
Regions     <- Regions[Regions$chrom == refChrom,]
print(head(Regions))

outFile <- paste0("polycomb_domains_counts.tab")
scoresAll <- data.frame()
head(scoresAll)

totalContacts <- read.table("obsContacts_per_sample.tab", header=T)
print(totalContacts)

outfile <- paste0("obsContacts_within_PcG_domains_vs_",refSample,"_in_",refChrom,".tab")

df <- read.table(outfile,header=T)
print(head(df))

summ <- df              %>%
    group_by(sample)  %>%
    reframe(n = n(), score = -0.1)

levels <- samples
my_comparisons <- list(c("LE_WT","larvae_DWT"),c("LE_WT","pupae_DWT"))

df$sample<- factor(df$sample, levels = levels)

statAnalysis <- compare_means(Log2FC_ContactsVsTotal_LE_WT ~ sample, data = df, method = "wilcox.test")
write.table(statAnalysis, file = "statAnalysis_normObsContacts_between_PcG_domains.tab", sep="\t", row.names=FALSE, quote=FALSE)

p <- ggplot(df, aes(x=sample, y=-Log2FC_ContactsVsTotal_LE_WT, fill=sample, alpha=0.2)) +
     geom_violin(position = position_dodge(1), trim=T, scale = "width") +
     geom_boxplot(position = position_dodge(1), width=0.1, fill="white", color="black", outlier.shape = NA) +
     theme(panel.background = element_rect(fill = NA),
           panel.grid.major.x = element_blank(),
           panel.grid.major.y = element_blank(),
           #axis.title = element_text(face="bold",size=24),
           axis.title = element_text(size=12),	   
           axis.text = element_text(face="bold",size=11),
           axis.text.x = element_text(angle=60, hjust=1),
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none") +
     labs(x="",y="Log2((LE_WT_intraTAD/LE_WT_validPairs)/(Mut_intraTAD/Mut_validPairs))",fill="") +     
     scale_color_manual(values=colors) +
     scale_fill_manual(values=colors) +
     scale_y_continuous(breaks=c(-0.5,-0.25,0,0.25)) +
     stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p", vjust = 0.1) +
     geom_point(data=df[df$domain1 == "dac",]) +     
     guides(fill="none")

pdf(paste0("obsContacts_within_PcG_domains_vs_",refSample,"_Fig2b.pdf"))
print(p)
dev.off()
