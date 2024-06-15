require(misha)
require(shaman)

require(doParallel)

### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source("./scripts/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

tssCoordinates <- gintervals.load("intervals.ucscCanTSS")
rownames(tssCoordinates) <- tssCoordinates$geneName

insTracks <- list(
	larvae_DWT = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w100kb_2kb',
	LD_gypsy1  = 'insulation.INS_hic_LD_gypsy1_merge_dm6_BeS_w100kb_2kb',
	LD_gypsy2  = 'insulation.INS_hic_LD_gypsy2Enh_merge_dm6_BeS_w100kb_2kb',
	LD_gypsy3  = 'insulation.INS_hic_LD_gypsy3_merge_dm6_SD_w100kb_2kb'
	)

refTrack <- "larvae_DWT"

chr <- 'chr2L'

cols <- rainbow(length(names(insTracks)))

insIntrv <- expand1D(tssCoordinates['dac',1:3],0.25e6)

### Insulation quantification
insQuantData <- list()
# 0. Compute insulation tracks using scripts/sandrine_computeInsulationTracks.R
# Retrieve the names of all insulation tracks 
allTracks <- gtrack.ls('insulation')
for(set in names(insTracks))
{
    tracks     <- allTracks[grep(set,allTracks)]
    trackNames <- gsub('_hic','',gsub('_merge_dm6_SD','',gsub('_merge_dm6_BeS','',gsub('insulation.','',tracks))))
    #print(set)
    #print(tracks)

    data <- gextract(tracks,insIntrv,iterator=2e3,colnames=trackNames)
    data[is.na(data)] <- 0
		
    for(t in trackNames)
    {
	data[,t] <- clipQuantile(-data[,t],0.999)
	data[,t] <- scaleData(data[,t],0,1,1e-9)		
    }
    insQuantData[[set]] <- data
}

pre1   = gintervals("chr2L",16422435,16423525)
pre2   = gintervals("chr2L",16485929,16486572)
gypsy1 = gintervals("chr2L",16461008,16461009)
gypsy2   =  gintervals("chr2L",16466059,16466060)
wingEnh  <- gintervals('chr2L',16465517,16465518)
gypsy3 = gintervals("chr2L",16478142,16478143)
boundary  <- gintervals('chr2L',16494000,16494001)
boundary1 <- gintervals('chr2L',16354000,16354001)
boundary2 <- gintervals('chr2L',16494000,16494001)
r0        <- gintervals("chr2L",16348000,16348001)
r1        <- gintervals("chr2L",16354000,16354001)
r2        <- gintervals("chr2L",16360000,16360001)
r3        <- gintervals("chr2L",16366000,16366001)
r4        <- gintervals("chr2L",16372000,16372001)
r5        <- gintervals("chr2L",16378000,16378001)
r6        <- gintervals("chr2L",16384000,16384001)
r7        <- gintervals("chr2L",16390000,16390001)
r8        <- gintervals("chr2L",16396000,16396001)
r9        <- gintervals("chr2L",16402000,16402001)
r10       <- gintervals("chr2L",16408000,16408001)
r11       <- gintervals("chr2L",16414000,16414001)
r12       <- gintervals("chr2L",16420000,16420001)

pre      <- gintervals('chr2L',16422514,16422515)
targetPoints <- list(gypsy1=gypsy1,gypsy2=gypsy2,gypsy3=gypsy3,r0=r0, r1=r1, r2=r2, r3=r3, r4=r4, r5=r5, r6=r6, r7=r7, r8=r8, r9=r9, r10=r10, r11=r11, r12=r12)
 
pdf('insulation_quantifications_Fig4d.pdf',width=12,height=9)

#print("")

par(mfrow=c(2,2),mar=c(8,5,1.5,2))
#print(names(insQuantData))
for(target in names(targetPoints))
{
    targetPos <- targetPoints[[target]]
    print(paste0("Analysing ",target," at ",targetPos$chrom,":",targetPos$start,"-",targetPos$end))
    maxneighbors = as.integer((targetPos$end - targetPos$start) / 2e3) + 1

    wtData <- gintervals.neighbors(targetPos,insQuantData[[refTrack]],maxneighbors=maxneighbors)[,grep('kb',colnames(insQuantData[[refTrack]]))+ncol(targetPos)]    

    npoints <- nrow(wtData)*ncol(wtData)
    #print(paste0("Number of points ",nrow(wtData)*ncol(wtData)))

    pvs <- c()
    pvalues <- c()

    for(set in names(insQuantData))
    {	

	mutData <- gintervals.neighbors(targetPos,insQuantData[[set]],maxneighbors=maxneighbors)[,grep('kb',colnames(insQuantData[[set]]))+ncol(targetPos)]
	pv <- -log(t.test(wtData,mutData)$p.value)
	pvalue <- t.test(wtData,mutData)$p.value

	pvs <- append(pvs,pv)
	pvalues <- append(pvalues,pvalue)
    }

    names(pvs) <- names(insQuantData)
    names(pvalues) <- names(insQuantData)

    barplot(pvs, main=paste0(target," ",npoints),las=2, ylim=c(0,12),col=cols, ylab='-log(Pval)')
    text('P=0.01',x=1,y=-log(0.01),pos=3)
    abline(h=-log(0.01))
    print(pvalues)
}
dev.off() 
quit()