require(reshape2) # cat tabulate to matrix
library(misha)
library(shaman)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
source('./scripts/auxFunctions.R')
options("scipen"=999)
colors <- rainbow(10)

mDBloc <- './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

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

totalContacts <- read.table("./scripts/obsContacts_per_sample.tab", header=T)
print(totalContacts)

outfile <- paste0("obsContacts_within_PcG_domains.tab")
if(!file.exists(outfile))
{
    for(sample in samples)
    {
        print(paste0("Analysing ",sample))

        result  <- data.frame()
        results <- data.frame()

        total <- totalContacts[totalContacts$sample == sample,]$count
        print(total)

        for(g in 1:nrow(Regions))
        {
            region <- Regions[g,]    

	    chrom1      <- region$chrom
    	    start1      <- region$start
    	    end1        <- region$end
	    domain1     <- region$domain
    	    print(paste(g,"of",nrow(Regions),chrom1,start1,end1,domain1,sep=" "), quote=FALSE)

	    for(d in 1:g)
	    {
	        domain <- Domains[d,]	    
	    
	        chrom2    <- domain$chrom
	        start2    <- domain$start
	        end2      <- domain$end
	        domain2   <- domain$domain
	        if(domain1 != domain2)
	        {
	            next
	        }
	        print(paste(d,"of",g,chrom2,start2,end2,domain2,sep=" "), quote=FALSE)	    

                ##### To plot scores #####
   	        interval <- gintervals.2d(chrom1,start1,end1,chrom2,start2,end2)
	        track <- obsData[[sample]]
	        print(track)
	        for(j in 1:length(track))
                {
                    vTrack <- gvtrack.create(paste0("obs",j),track[j],"area")
                }
	        data <- gextract(paste0("obs",1:length(track)),interval,iterator=interval)

	        contacts <- sum(data[,paste0("obs",1:length(track))])
	    
	        result <- data.frame(sample,chrom1,start1,end1,domain1,chrom2,start2,end2,domain2,contacts,contacts/total,total)
	        if(nrow(results) == 0)
	        {
	            results <- result
	        } else {
	            results <- rbind(results,result)
	        }
	    
	    } # Close cycle over d
        } # Close cycle over r

        colnames(results) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","contacts","contacts_vs_total","total")
        write.table(as.matrix(results),file=outfile,sep="\t",quote=F,row.names=F,col.names=F,append=T)

    } # Close cycle over sample
}

df_obsContacts <- read.table(outfile,header=F)
colnames(df_obsContacts) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","contacts","contacts_vs_total","total")
print(head(df_obsContacts))

result <- data.frame()
results <- data.frame()

outfile <- paste0("obsContacts_within_PcG_domains_vs_",refSample,"_in_",refChrom,".tab")
if(!file.exists(outfile))
{
    for(domain in unique(df_obsContacts$domain1))
    {
        print(domain)

        domainContacts <- df_obsContacts[df_obsContacts$domain1 == domain,]
    
        chrom1     <- unique(domainContacts$chrom1)
        start1     <- unique(domainContacts$start1)
        end1       <- unique(domainContacts$end1)
        chrom2     <- unique(domainContacts$chrom2)
        start2     <- unique(domainContacts$start2)
        end2       <- unique(domainContacts$end2)    

        refContacts <- domainContacts[domainContacts$sample == refSample,]$contacts
        refRatio    <- domainContacts[domainContacts$sample == refSample,]$contacts_vs_total
 
        for(sample in unique(domainContacts$sample))
        {   
 
            obsContacts <- domainContacts[domainContacts$sample == sample,]$contacts
            obsRatio    <- domainContacts[domainContacts$sample == sample,]$contacts_vs_total
	    result <- data.frame(sample,chrom1,start1,end1,domain,chrom2,start2,end2,domain,refContacts,refRatio,obsContacts,obsRatio,log(obsContacts/refContacts)/log(2),log(obsRatio/refRatio)/log(2))

	    if(nrow(results) == 0)
            {
	        results <- result
	    } else {
	        results <- rbind(results,result)
	    }
	
        }
    }

    colnames(results) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","refContacts","refRatio","obsContacts","obsRatio",paste0("Log2FC_Contacts_",refSample),paste0("Log2FC_ContactsVsTotal_",refSample))
    write.table(as.matrix(results),file=outfile,sep="\t",quote=F,row.names=F)
}
