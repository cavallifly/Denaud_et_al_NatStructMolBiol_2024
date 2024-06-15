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

diffData <- list(
	LD_gypsy1_k250_kexp250 = "hic.diffMaps.hic_larvae_DWT_merge_dm6_BeS_vs_hic_LD_gypsy1_merge_dm6_BeS_chr2L_score_k250_kexp250_1Mb",	
	LD_gypsy2_k250_kexp250 = "hic.diffMaps.hic_larvae_DWT_merge_dm6_BeS_vs_hic_LD_gypsy2Enh_merge_dm6_BeS_chr2L_score_k250_kexp250_1Mb",		
	LD_gypsy3_k250_kexp250 = "hic.diffMaps.hic_larvae_DWT_merge_dm6_BeS_vs_hic_LD_gypsy3_merge_dm6_SD_chr2L_score_k250_kexp250_1Mb"	
)
print(diffData)

resolution <- c(3e3)

chr <- "chr2L"

# Genomic locations of interest
# From Bernd
#Gypsy 1 2L : 16461008
#Gypsy 2 2L : 16466059
#Gypsy 3 2L : 16478142
#PRE1 : 2L : 16422435..16423525
#PRE2 : 2L : 16485929..16486572

pre    = gintervals("chr2L",16422500,16422501)
pre1   = gintervals("chr2L",16422435,16423525)
pre2   = gintervals("chr2L",16485929,16486572)
gypsy1 = gintervals("chr2L",16461008,16461009)
gypsy2 = gintervals("chr2L",16466059,16466060)
gypsy3 = gintervals("chr2L",16478142,16478143)

preAnnotation <- g2d(rbind(expand1D(pre1,1e3),expand1D(pre2,1e3))[,1:3])
gypsy1Annotation <- g2d(rbind(expand1D(gypsy1,1e3))[,1:3])
gypsy2Annotation <- g2d(rbind(expand1D(gypsy2,1e3))[,1:3])
gypsy3Annotation <- g2d(rbind(expand1D(gypsy3,1e3))[,1:3])

tags <- c("k250_kexp250")

for(annotated in c("FALSE"))
{
    for(tag in tags)
    {
        print(tag)
	res <- resolution	

	trackData <- diffData[grep(tag,diffData)]
	print(trackData)

        start <- 16354000
	end   <- 16500000	    

        intrv <- gintervals(chr,start,end)
        exIntrv <- giterator.intervals(intervals=g2d(intrv),iterator=c(res,res))

        nr <- 1
	nc <- length(names(trackData))
	slotSize <- 800

	diffCol <- colorRampPalette(c("darkred","lightgrey","darkblue"))(200)

	png(file=paste0("diffMaps_",res,"bp_",tag,"_Fig4b_annot_",annotated,".png"),width=nc*slotSize,height=nr*slotSize)	
        layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
        par(mai=rep(0.45,4))

        for(set in names(trackData))
        {

            track <- trackData[[set]]
	    vTrack <- gvtrack.create("score",track,"avg")
	    data <- gextract("score",exIntrv,iterator=exIntrv,colnames="score")
	    data[is.na(data)] <- 0
	    outMatrix <- paste0("diffMaps_",res,"bp_",tag,"_Fig4b_",set,".tsv")
	    write.table(data[,c(1,2,3,4,5,6,7)], file = outMatrix, sep="\t", row.names=FALSE, quote=FALSE)

	    mtx <- acast(data,start1~start2,value.var="score")

	    counts <- gextract(track,g2d(intrv),colnames="score")
	    regReads <- nrow(counts)/2
	    print(set)

	    setName = gsub(tag,'',gsub('hic_','',gsub('_score','',gsub('_merge_dm6_BeS','',gsub('hic.diffMaps.hic_','',set)))))
	    print(setName)

	    title <- paste0(setName," - ", signif(regReads/1000,3),"k in region (",signif(min(data$score),2),":",signif(max(data$score),2),")")
	    image(mtx, main=title, useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=diffCol)
	    axis(1,at=c(0,1.5,1),labels=c(signif(start/1e6,3),chr,signif(end/1e6,3)),cex.axis=2.6)

	    if(annotated == "TRUE")
	    {
	        # PREs annotation on the diagonal
	        plotRectangles2D(preAnnotation,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	        # PRE loop annotation
	        PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
	        plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)

	        # gypsy insertion annotation (sample specific)
	        if(setName == "LD_gypsy1")
	        {
	            plotRectangles2D(gypsy1Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	        if(setName == "LD_gypsy2")
	        {
	            plotRectangles2D(gypsy2Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	        if(setName == "LD_gypsy3")
	        {
	            plotRectangles2D(gypsy3Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	    }
	}
	dev.off()

	png(file=paste0("diffMaps_",res,"bp_",tag,"_Fig4b_annot_",annotated,"_toModify.png"),width=nc*slotSize,height=nr*slotSize)	
	layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
	par(mai=rep(0.45,4))

	plotData <- list()
	    
	for(set in names(trackData))
        {

	    track <- trackData[[set]]
	    vTrack <- gvtrack.create("score",track,"avg")
	    data <- gextract("score",exIntrv,iterator=exIntrv,colnames="score")
	    data[is.na(data)] <- 0
	    mtx <- acast(data,start1~start2,value.var="score")

            counts <- gextract(track,g2d(intrv),colnames="score")
	    regReads <- nrow(counts)/2
	    print(set)

	    setName = gsub(tag,'',gsub('hic_','',gsub('_score','',gsub('_merge_dm6_BeS','',gsub('hic.diffMaps.hic_','',set)))))
	    print(setName)

	    title <- paste0(setName," - ", signif(regReads/1000,3),"k in region (",signif(min(data$score),2),":",signif(max(data$score),2),")")
	    image(mtx, main="", useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=diffCol) #, zlim=c(-75,75))		

	    if(annotated == "TRUE")
	    {
	        # PREs annotation on the diagonal
	        plotRectangles2D(preAnnotation,g2d(intrv),lwd=6,image=TRUE,col='yellow')
	        # PRE loop annotation
	        PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
	        plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)

	        # gypsy insertion annotation (sample specific)
	        if(setName == "LD_gypsy1")
	        {
	            plotRectangles2D(gypsy1Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	        if(setName == "LD_gypsy2")
	        {
	            plotRectangles2D(gypsy2Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	        if(setName == "LD_gypsy3")
	        {
	            plotRectangles2D(gypsy3Annotation,g2d(intrv),lwd=6,image=TRUE,col='cyan')
	        }
	    }
	}
	dev.off()

	# Plot the colorbar once for all
	pdf(paste0('diffMaps_colorbar.pdf'),width=2,height=8)
	labels <- c("-100","-50","0","50","100")
	x <- length(diffCol)
	plot(y=1:x,x=rep(0,x),col=diffCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
	sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
	dev.off()
    }
}