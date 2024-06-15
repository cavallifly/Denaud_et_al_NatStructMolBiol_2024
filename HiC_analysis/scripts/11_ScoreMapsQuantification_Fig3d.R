library(misha)
library(shaman)
library(doParallel)

### misha working DB
mDBloc <-  "./mishaDB/trackdb/"
db <- "dm6"
dbDir <- paste0(mDBloc,db,"/")
gdb.init(dbDir)
gdb.reload()

source("./scripts/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

k    <- 250
kexp <- 2*k
pTarget <- "larvae_DWT"
scoreData <- list(
	larvae_DWT = paste0("hic.hicScores_larvae_DWT_merge_dm6_BeS.hicScores_larvae_DWT_merge_dm6_BeS_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb"),
	LD_DRE     = paste0("hic.hicScores_LD_DRE_merge_dm6_BeS.hicScores_LD_DRE_merge_dm6_BeS_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb")
        )

scoreTracks <- scoreData
print(scoreTracks)

# Genomic locations of interest
pre      = gintervals("chr2L",16422500,16422501)
pre1     = gintervals("chr2L",16422514,16422515)
pre2     = gintervals("chr2L",16485929,16485930)
wingEnh  = gintervals("chr2L",16465517,16465518)
gypsy    = gintervals("chr2L",16460991,16461011)
gypsyEnh = gintervals("chr2L",16466064,16466065)
boundary = gintervals("chr2L",16494000,16494001)
boundary1 <- gintervals('chr2L',16340000,16340001) ###
boundary2 <- gintervals('chr2L',16494000,16494001)
viewpoints <- list(
	   pre1  = gintervals("chr2L",16422514,16422515),
	   pre2  = gintervals("chr2L",16485929,16485930)
	   )
vPoints <- do.call(rbind.data.frame, viewpoints)

####################### Quantify Significance
vPairs <- build2D(vPoints,vPoints,diag=FALSE, minDist=10e3)

scoreTracks <- scoreData
pdf(paste0("PRE_loop_quantifications_k",k,"_kexp",kexp,"_Fig3d_comp_",pTarget,".pdf"),width=14)

layout(matrix(1:3,ncol=3,byrow=TRUE),respect=FALSE)

res <- 3.0e3
testedIntrvs <- c()
pointData <- list()
for(i in 1:nrow(vPoints))
{
    for(j in 1:nrow(vPoints))
    {
        if(i == j){next}
    	if(vPoints[i,"start"] - vPoints[j,"start"] < 10e3){next}

    	intrv <- gintervals.2d(vPoints[i,"chrom"],vPoints[i,"start"],vPoints[i,"end"],vPoints[j,"chrom"],vPoints[j,"start"],vPoints[j,"end"])
    	pointName <- paste(rownames(vPoints)[i],rownames(vPoints)[j],sep="_")

    	testedIntrvs <- rbind(testedIntrvs,intrv)

    	print(paste(rownames(vPoints)[i],rownames(vPoints)[j]))
    	print(pointName)
    	print(intrv)
    	intrv <- expand2D(intrv,res)

    	for(set in names(scoreTracks))
    	{
	    setName = gsub('hic.','',gsub('_merge_dm6_SD','',gsub('_merge_dm6_BeS','',gsub('hicScores_','',set))))
            print(setName)
	    data <- gextract(scoreTracks[[set]],intrv,colnames="score")
	    pointData[[pointName]][[setName]] <- data$score
	    outMatrix <- paste0("scoreMapsk",k,"kexp",kexp,"_r",res,"bp_Fig3d_",set,".tsv")
	    write.table(data[,c(1,2,3,4,5,6,7)], file = outMatrix, sep="\t", row.names=FALSE, quote=FALSE)
    	}
    	jitterPlot(pointData[[pointName]],pTarget=pTarget,main=pointName, ylim=c(0,100),jitWidth=.75,cols=c("darkslateblue",rep("darkslategrey",length(names(scoreTracks))-1)))
    	####
    }
}
dev.off()
