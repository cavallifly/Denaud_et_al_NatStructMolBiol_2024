library(misha)
library(shaman)
library(doParallel)

##########################################
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
zmin <- -100
zmax <-  100
scoreData <- list(
        larvae_DWT = paste0("hic.hicScores_larvae_DWT_merge_dm6_BeS.hicScores_larvae_DWT_merge_dm6_BeS_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb"),
        LD_gypsy1  = paste0("hic.hicScores_LD_gypsy1_merge_dm6_BeS.hicScores_LD_gypsy1_merge_dm6_BeS_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb"),
        LD_gypsy2  = paste0("hic.hicScores_LD_gypsy2Enh_merge_dm6_BeS.hicScores_LD_gypsy2Enh_merge_dm6_BeS_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb"),
        LD_gypsy3  = paste0("hic.hicScores_LD_gypsy3_merge_dm6_SD.hicScores_LD_gypsy3_merge_dm6_SD_chr2L_10_20Mb_score_k",k,"_kexp",kexp,"_step1Mb")
        )

resolution <- c(3.0e3)

chr <- "chr2L"

CG5888prom = gintervals("chr2L",16448468,16448469)
RE         = gintervals("chr2L",16465517,16465518)
pre1   = gintervals("chr2L",16422435,16423525)
pre2   = gintervals("chr2L",16485929,16486572)
gypsy1 = gintervals("chr2L",16461008,16461009)
gypsy2 = gintervals("chr2L",16466059,16466060)
gypsy3 = gintervals("chr2L",16478142,16478143)

CG5888Annotation <- g2d(rbind(expand1D(CG5888prom,1e3))[,1:3])
REAnnotation     <- g2d(rbind(expand1D(RE,1e3))[,1:3])
preAnnotation    <- g2d(rbind(expand1D(pre1,1e3),expand1D(pre2,1e3))[,1:3])
gypsy1Annotation <- g2d(rbind(expand1D(gypsy1,1e3))[,1:3])
gypsy2Annotation <- g2d(rbind(expand1D(gypsy2,1e3))[,1:3])
gypsy3Annotation <- g2d(rbind(expand1D(gypsy3,1e3))[,1:3])

for(annotated in c("FALSE"))
{

    res <- resolution

    start <- 16321000
    end   <- 16512000

    intrv <- gintervals(chr,start,end)
    exIntrv <- giterator.intervals(intervals=g2d(intrv),iterator=c(res,res))

    ### Plot Score maps
    nr <- 1
    slotSize <- 800

    nc <- length(names(scoreData))

    png(file=paste0("scoreMapsk",k,"kexp",kexp,"_r",res,"bp_Fig4a_annot_",annotated,".png"),width=nc*slotSize,height=nr*slotSize)
    print(paste0("Plotting map between ",start,"-",end))
    
    scoreCol <- bonevScale()

    layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
    par(mai=rep(0.45,4))

    plotData <- list()

    for(set in names(scoreData))
    {
        track  <- scoreData[[set]]
	data   <- gextract(track,exIntrv,iterator=exIntrv,colnames="score")
	data[is.na(data)] <- 0
	outMatrix <- paste0("scoreMapsk",k,"kexp",kexp,"_r",res,"bp_Fig4a_",set,".tsv")
	write.table(data[,c(1,2,3,4,5,6,7)], file = outMatrix, sep="\t", row.names=FALSE, quote=FALSE)

	mtx    <- acast(data,start1~start2,value.var="score")

	counts <- gextract(track,g2d(intrv),colnames="score")
	regReads <- nrow(counts)/2

	title <- paste0(set," - ", signif(regReads/1000,3),"k in region (",signif(min(data$score),2),":",signif(max(data$score),2),")")

	setName = gsub('hic_','',gsub('_score','',gsub('_merge_dm6_BeS','',gsub('hic.diffMaps.hic_','',set))))

	image(mtx, main=title, useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=scoreCol, zlim=c(zmin,zmax))
	axis(1,at=c(0,0.5,1),labels=c(signif(start/1e6,3),chr,signif(end/1e6,3)),cex.axis=2.6)
	if (annotated == "TRUE")
	{
            # Annotation on the diagonal
	    plotRectangles2D(preAnnotation   ,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	    plotRectangles2D(CG5888Annotation,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	    plotRectangles2D(REAnnotation    ,g2d(intrv),lwd=6,image=TRUE,col='yellow')	

	    # PRE loop annotation
	    PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
	    plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)

	    # gypsy insertion annotation (sample specific)
	    if(setName == "LD_gypsy1")
	    {
		plotRectangles2D(gypsy1Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
	    if(setName == "LD_gypsy2")
	    {
		plotRectangles2D(gypsy2Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
	    if(setName == "LD_gypsy3")
	    {
                plotRectangles2D(gypsy3Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
 	}
    }
    dev.off()

    png(file=paste0("scoreMapsk",k,"kexp",kexp,"_r",res,"bp_Fig4a_annot_",annotated,"_toModify.png"),width=nc*slotSize,height=nr*slotSize)
    print(paste0("Plotting map between ",start,"-",end))

    scoreCol <- bonevScale()

    layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
    par(mai=rep(0.45,4))

    for(set in names(scoreData))
    {
        track  <- scoreData[[set]]
	data   <- gextract(track,exIntrv,iterator=exIntrv,colnames="score")
	data[is.na(data)] <- 0
	mtx    <- acast(data,start1~start2,value.var="score")

	counts <- gextract(track,g2d(intrv),colnames="score")
	regReads <- nrow(counts)/2

	setName = gsub('hic_','',gsub('_score','',gsub('_merge_dm6_BeS','',gsub('hic.diffMaps.hic_','',set))))
	image(mtx, main="", useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=scoreCol, zlim=c(zmin,zmax))
	if (annotated == "TRUE")
	{
            # Annotation on the diagonal
	    plotRectangles2D(preAnnotation   ,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	    plotRectangles2D(CG5888Annotation,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	    plotRectangles2D(REAnnotation    ,g2d(intrv),lwd=6,image=TRUE,col='yellow')	

	    # PRE loop annotation
	    PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
	    plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)

	    # gypsy insertion annotation (sample specific)
	    if(setName == "LD_gypsy1")
	    {
		plotRectangles2D(gypsy1Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
	    if(setName == "LD_gypsy2")
	    {
                plotRectangles2D(gypsy2Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
	    if(setName == "LD_gypsy3")
	    {
		plotRectangles2D(gypsy3Annotation,g2d(intrv),lwd=6,image=TRUE,col='darkgreen')
	    }
 	}
    }
    dev.off()

    # Plot once for all the colorbar
    pdf(paste0('scoreMaps_colorbar.pdf'),width=2,height=8)
    labels <- c("-100","-50","0","50","100")
    x <- length(scoreCol)
    plot(y=1:x,x=rep(0,x),col=scoreCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
    sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
    dev.off()

}