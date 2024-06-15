library(misha)
library(shaman)

source("./scripts/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

### misha working DB
mDBloc <-  "./mishaDB/trackdb/"
db     <- "dm6"
dbDir  <- paste0(mDBloc,db,"/")
gdb.init(dbDir)
gdb.reload()

tssCoordinates           <- gintervals.load("intervals.ucscCanTSS")
rownames(tssCoordinates) <- tssCoordinates$geneName

insTracks <- list(
        larvae_DWT_100  = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w100kb_r2kb',
        larvae_DWT_150  = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w150kb_r2kb',
        larvae_DWT_200  = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w200kb_r2kb',
        larvae_DWT_250  = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w250kb_r2kb',
        larvae_DWT_300  = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w300kb_r2kb',	
	LD_DRE_100      = 'insulation.INS_hic_LD_DRE_merge_dm6_BeS_w100kb_r2kb',
	LD_DRE_150      = 'insulation.INS_hic_LD_DRE_merge_dm6_BeS_w150kb_r2kb',
	LD_DRE_200      = 'insulation.INS_hic_LD_DRE_merge_dm6_BeS_w200kb_r2kb',
	LD_DRE_250      = 'insulation.INS_hic_LD_DRE_merge_dm6_BeS_w250kb_r2kb',
	LD_DRE_300      = 'insulation.INS_hic_LD_DRE_merge_dm6_BeS_w300kb_r2kb'
        )
refSet = "larvae_DWT"

# Genomic locations of interest
pre1     = gintervals("chr2L",16422514,16422515)
gypsy1   = gintervals("chr2L",16461008,16461011)
gypsy2   = gintervals("chr2L",16466059,16466060)
gypsy3   = gintervals("chr2L",16478142,16478143)
boundary1 <- gintervals("chr2L",16354000,16354001)
boundary2 <- gintervals("chr2L",16494000,16494001)
dacTSS   = tssCoordinates["dac",1:3]
wingEnh  <- gintervals('chr2L',16465517,16465518)

r0        <- gintervals("chr2L",16348000,16348001) ; print(paste0("R0 ",r0))
r1        <- gintervals("chr2L",16354000,16354001) ; print(paste0("R1 ",r0))
r2        <- gintervals("chr2L",16360000,16360001) ; print(paste0("R2 ",r0))
r3        <- gintervals("chr2L",16366000,16366001) ; print(paste0("R3 ",r0))
r4        <- gintervals("chr2L",16372000,16372001) ; print(paste0("R4 ",r0))
r5        <- gintervals("chr2L",16378000,16378001) ; print(paste0("R5 ",r0))
r6        <- gintervals("chr2L",16384000,16384001) ; print(paste0("R6 ",r0))
r7        <- gintervals("chr2L",16390000,16390001) ; print(paste0("R7 ",r0))
r8        <- gintervals("chr2L",16396000,16396001) ; print(paste0("R8 ",r0))
r9        <- gintervals("chr2L",16402000,16402001) ; print(paste0("R9 ",r0))
r10       <- gintervals("chr2L",16408000,16408001) ; print(paste0("R10",r0))
r11       <- gintervals("chr2L",16414000,16414001) ; print(paste0("R11",r0))
r12       <- gintervals("chr2L",16420000,16420001) ; print(paste0("R12",r0))

chr <- "chr2L"
center <- 16485998
insIntrv <- gintervals(chr, center - 165000, center + 25000)

# Retrieve the insulation tracks from mishaDB
print("Get the min, mean, stddev, and max of the insulation tracks")
insData <- list()
for(set in names(insTracks))
{
    print(set, row.names=F, quote=F)
    track <- insTracks[[set]]
    data <- gextract(track,insIntrv,iterator=3e3,colnames="INS")	
    data$INS  <- clipQuantile(-data$INS,0.999)
    data$INS  <- scaleData(data$INS,0,1,1e-9)
    
    trackName  <- gsub("_hic","",gsub("_merge_dm6_SD","",gsub("_merge_dm6_BeS","",gsub("insulation.","",set))))
    trackName1 <- unlist(strsplit(trackName,"_"))[1:2]
    trackName1 <- paste(trackName1,collapse="_")
    window     <- paste0("w",unlist(strsplit(trackName,"_"))[3])
    if(length(insData[[trackName1]]) == 0)
    {
        insData[[trackName1]] <- data[c(1,2,3,4)]
	colnames(insData[[trackName1]]) <- c("chrom","start","end",window)
    } else {
        cNames <- c(colnames(insData[[trackName1]]),window)
        insData[[trackName1]] <- cbind(insData[[trackName1]],data$INS)
	colnames(insData[[trackName1]]) <- cNames
    }

}

print("")
for(trackName in names(insData))
{
    print(trackName)
    print(head(insData[[trackName]]))

    insData[[trackName]]$min    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,min)    
    insData[[trackName]]$mean   <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,mean)
    insData[[trackName]]$stddev <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,sd)
    insData[[trackName]]$max    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,max)

    outProfile <- paste0("insulationProfiles_Fig3e_",trackName,".tsv")
    write.table(insData[[trackName]], file = outProfile, sep="\t", row.names=FALSE, quote=FALSE)

    print("", row.names=F, quote=F)
}
#print(head(insData))


# Obtain the location on the plot of the GLoI
r0StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r0$start))/(insIntrv$end - insIntrv$start))
r1StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r1$start))/(insIntrv$end - insIntrv$start))
r2StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r2$start))/(insIntrv$end - insIntrv$start))
r3StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r3$start))/(insIntrv$end - insIntrv$start))
r4StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r4$start))/(insIntrv$end - insIntrv$start))
r5StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r5$start))/(insIntrv$end - insIntrv$start))
r6StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r6$start))/(insIntrv$end - insIntrv$start))
r7StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r7$start))/(insIntrv$end - insIntrv$start))
r8StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r8$start))/(insIntrv$end - insIntrv$start))
r9StartPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r9$start))/(insIntrv$end - insIntrv$start))
r10StartPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r10$start))/(insIntrv$end - insIntrv$start))
r11StartPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r11$start))/(insIntrv$end - insIntrv$start))
r12StartPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r12$start))/(insIntrv$end - insIntrv$start))
gypsy1StartPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy1$start))/(insIntrv$end - insIntrv$start))
gypsy2StartPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy2$start))/(insIntrv$end - insIntrv$start))
gypsy3StartPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy3$start))/(insIntrv$end - insIntrv$start))
r0EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r0$end))/(insIntrv$end - insIntrv$start))
r1EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r1$end))/(insIntrv$end - insIntrv$start))
r2EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r2$end))/(insIntrv$end - insIntrv$start))
r3EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r3$end))/(insIntrv$end - insIntrv$start))
r4EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r4$end))/(insIntrv$end - insIntrv$start))
r5EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r5$end))/(insIntrv$end - insIntrv$start))
r6EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r6$end))/(insIntrv$end - insIntrv$start))
r7EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r7$end))/(insIntrv$end - insIntrv$start))
r8EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r8$end))/(insIntrv$end - insIntrv$start))
r9EndPos     <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end -  r9$end))/(insIntrv$end - insIntrv$start))
r10EndPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r10$end))/(insIntrv$end - insIntrv$start))
r11EndPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r11$end))/(insIntrv$end - insIntrv$start))
r12EndPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - r12$end))/(insIntrv$end - insIntrv$start))
gypsy1EndPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy1$end))/(insIntrv$end - insIntrv$start))
gypsy2EndPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy2$end))/(insIntrv$end - insIntrv$start))
gypsy3EndPos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy3$end))/(insIntrv$end - insIntrv$start))
dacPos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - dacTSS$start))/(insIntrv$end - insIntrv$start))
pre1Pos    <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - pre1$start))/(insIntrv$end - insIntrv$start))
wingPos   <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - wingEnh$start))/(insIntrv$end - insIntrv$start))
gypsy3Pos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - gypsy3$start))/(insIntrv$end - insIntrv$start))
boundary1Pos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - boundary1$start))/(insIntrv$end - insIntrv$start))
boundary2Pos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - boundary2$start))/(insIntrv$end - insIntrv$start))

cols     <- rainbow(length(names(insMean)))
fillCols <- rainbow(length(names(insMean)))
cols[[1]]     <- "black"
fillCols[[1]] <- rgb(0,     0, 0, max = 255, alpha = 50, names = "black50")
cols[[2]]     <- "#FFc000"
fillCols[[2]] <- rgb(255, 192, 0, max = 255, alpha = 50, names = "orange50")
print(cols)

print(cols) ; print(fillCols)

for(set in names(insData))
{

    if(set == "larvae_DWT")
    {
        next
    }

    if(set == "LD_DRE")
    {
        pdf(paste0("insulationProfiles_Fig3e.pdf"),width=19,height=7)    
	plot(insData[[1]]$mean, type="l", ylim=c(0,1),xlim=c(0,nrow(insData[[1]])+nrow(insData[[1]])*0.15),xlab="", ylab="INSULATION", xaxt="n", frame=F)
	polygon(c(1:nrow(insData[[refSet]]), nrow(insData[[refSet]]):1), c(insData[[refSet]]$mean-insData[[refSet]]$stddev,rev(insData[[refSet]]$mean+insData[[refSet]]$stddev)), col = fillCols[which(names(insData) == refSet)], border=F)

	axis(1,c(insIntrv$start,"PRE1","TAD1",insIntrv$end),at=c(1,pre1Pos,boundary1Pos,nrow(insData[[1]])))	
        legend("topright",legend=names(insData[c(1,2)]), col=cols[c(1,2)], pch=19)
    }


    print(set)

    layout(matrix(1:1,ncol=1,byrow=TRUE),respect=FALSE)
    par(mai=rep(1.5,4), cex=1.5)

    abline(v=dacPos, col=4, lwd=3, lty=2)
    abline(v=prePos, col=8, lwd=3, lty=2)
    abline(v=boundary1Pos, col=5, lwd=3, lty=2)
    abline(v=boundary2Pos, col=5, lwd=3, lty=2)
    abline(v=wingPos, col=2, lwd=3, lty=2)

    lines(insData[[set]]$mean, type="l", col=cols[which(names(insData) == set)], lwd=3)
    polygon(c(1:nrow(insData[[set]]), nrow(insData[[set]]):1), c(insData[[set]]$mean-insData[[set]]$stddev,rev(insData[[set]]$mean+insData[[set]]$stddev)), col = fillCols[which(names(insData) == set)], border=F)

    dev.off()    
}
