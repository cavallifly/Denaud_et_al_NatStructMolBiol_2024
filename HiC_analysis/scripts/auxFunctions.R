require(misha)
require(shaman)
require(plyr)
require(reshape2)
require(pheatmap)
require(plotrix)

### BEGIN : scaleData ###
scaleData <- function(data,from=0,to=1,round=0.01)
{
    xM <- max(data,na.rm=T)
    xm <- min(data,na.rm=T)
    data <- sapply(data,function(x){(x-xm)/(xM-xm)})
	
    data <- round_any((data*(to-from)),round)
    data <- data + from
	
    data[data > to] <- to
    data[data < from] <- from
	
    return(data)
}
### END : scaleData ###

### BEGIN : g2d ###
g2d <- function(intrv)
{
    return(gintervals.2d(intrv$chrom,intrv$start,intrv$end,intrv$chrom,intrv$start,intrv$end))
}
### BEGIN : g2d ###

### BEGIN : jitterPlot ###
jitterPlot <- function(data,cols=NA,main='',ylim=NA,pTarget=NA,points=TRUE,cexPV=1,valPV=TRUE,jitWidth=0.5,pvalues=TRUE,test='wilcox',alternative='two.sided')
{
    if(!is.list(data))
    {
        print('PLEASE PROVIDE A NAMED LIST')
	return()
    }
	
    options(scipen=0)
	
    bCols <- rep('grey40', length(data))
    if(!is.na(pTarget))
    {
        bCols[which(names(data) == pTarget)] <- 'blue'
    }
	
    if(is.na(ylim[1]))
    {
	plt <- boxplot(data, notch=T, frame=F, outline=F, main=main, col=bCols,las=2)
	yMax <- max(plt$stats)*0.85
    }else{
        boxplot(data, notch=T, frame=F, outline=F, main=main, col=bCols,ylim=ylim,las=2)
	yMax <- ylim[2]*0.85
    }

    abline(h=0, lwd=1.5,col='black',lty=2)
    for (i in 1:length(data))
    {
  	Y <- data[[i]]
	X <- rep(i,length(Y))
	ifelse(!is.na(cols),col <- cols[i],col <- i+1)
	if(points ==TRUE)
	{
	    points(jitter(X,length(data)/i*jitWidth),Y,pch=19,cex=0.25, col=col)
	    text(x=i+0.,y=3,length(Y),cex=cexPV)			
	}
		
	if(pvalues == TRUE)
	{
	    pPlot <- FALSE
	    if(is.na(pTarget))
	    {		
		if(i < length(data))
		{
		    if(test == 'wilcox')
		    {
		        pval <- wilcox.test(data[[i]],data[[i+1]],alternative=alternative)$p.value
			print(wilcox.test(data[[i]],data[[i+1]],alternative=alternative))
		    }else{
			pval <- ks.test(data[[i]],data[[i+1]],alternative=alternative)$p.value
		    }
		    pPlot <- TRUE
		}
	    }else{
	    	if(test == 'wilcox')
		{
		    pval <- wilcox.test(data[[i]],data[[pTarget]],alternative=alternative)$p.value
		    print(  wilcox.test(data[[i]],data[[pTarget]],alternative=alternative))
		}else{
		    pval <- ks.test(data[[i]],data[[pTarget]],alternative=alternative)$p.value
		}
		if(names(data)[i] != pTarget){pPlot <- TRUE}
		i <- i-0.5
	    }
		
	    if(pPlot == TRUE)
	    {
		if(valPV == TRUE)
		{
		    text(x=i+0.5,y=yMax*1.1,signif(pval,2),cex=cexPV)
		    ySig <- yMax
		}else{
		    ySig <- yMax*1.1
		}
			
		if(pval <= 1e-6)
		{
		    text(x=i+0.5,y=ySig,'***',cex=cexPV*1.1,col=2)
		}else if(pval <= 1e-3){
		    text(x=i+0.5,y=ySig,'**',cex=cexPV*1.1,col=4)
		}else if(pval <= 0.01){
		    text(x=i+0.5,y=ySig,'*',cex=cexPV*1.1)
		}else{
		    text(x=i+0.5,y=ySig,'NS',cex=cexPV*0.8)
		}
				
            }
	}
    }
}
### END : jitterPlot ###

### BEGIN : clipQuantile ###
clipQuantile <- function(set, thr, up=TRUE, down=TRUE)
{
    set <- as.numeric(as.character(set))

    upLim   <- quantile(set,thr,na.rm=TRUE)
    downLim <- quantile(set,1-thr,na.rm=TRUE)
	
    set[set <= downLim] <- downLim
    set[set >= upLim] <- upLim
	
    return(set)
}
### END : clipQuantile ###

### BEGIN : bonevScale ###
bonevScale <- function()
{
    palette.breaks = function(n, colors, breaks)
    {
	colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
	for(i in 2:(length(colors)))
	{
	    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
	}
	colspec = c(colspec,colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)]))
	colspec
    }
    col.pbreaks <- c(20,35,50,65,75,85,95)
    col.pos <- palette.breaks(100 , c("lightgrey","lavenderblush2","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
    col.nbreaks <- c(20,35,50,65,75,85,95)
    col.neg <- rev(palette.breaks(100 , c("lightgrey", "powderblue", "cornflowerblue", "blue","blueviolet", "mediumpurple4", "black"), col.nbreaks ))
    scoreCol <- c(col.neg, col.pos)
    return(scoreCol)
}
### END : bonevScale ###

### BEGIN : build2D ###
build2D <- function(set1,set2,names1='',names2='',minDist=0,maxDist=10e10,band=TRUE,diag=TRUE,retNames=TRUE,trans=FALSE)
{

    ####################################################

    chrs <- unique(c(as.character(set1$chrom),as.character(set2$chrom)))
    ints <- c()
    for(chr in chrs)
    {
	d1 <- set1[as.character(set1$chrom) == chr,]
	d2 <- set2[as.character(set2$chrom) == chr,]
		
	if(nrow(d1) == 0 | nrow(d2) == 0){next}
		
	d_inds     = as.data.frame(expand.grid(1:nrow(d1), 1:nrow(d2)))
	d_ints_for = data.frame(chrom1=d1[d_inds[,1],'chrom'], start1=d1[d_inds[,1],'start'], end1=d1[d_inds[,1],'end'], chrom2=d2[d_inds[,2],'chrom'], start2=d2[d_inds[,2],'start'], end2=d2[d_inds[,2],'end'])
	d_ints_rev = data.frame(chrom1=d2[d_inds[,2],'chrom'], start1=d2[d_inds[,2],'start'], end1=d2[d_inds[,2],'end'], chrom2=d1[d_inds[,1],'chrom'], start2=d1[d_inds[,1],'start'], end2=d1[d_inds[,1],'end'])
		
	d_ints <- rbind(d_ints_for,d_ints_rev)
	d_ints <- d_ints[!duplicated(paste(d_ints$start1,d_ints$end1,d_ints$start2,d_ints$end2,sep='_')),]
		
	##################
	d_ints <- gintervals.canonic(d_ints)
		
	if(is.null(ints))
	{
	    ints <- d_ints
	}else{
	    ints <- rbind(ints,d_ints)
	}
    }

    if(is.null(ints))
    {
        return(ints)
    }

    if(band==TRUE)
    {
        ints <- ints[ints$start2 >= ints$start1,]
    }
	
    ints$dist <- abs(ints$start2 - ints$start1)
    ints <- ints[ints$dist >= minDist & ints$dist <= maxDist,]
	
    if(diag==FALSE)
    {
        ints <- ints[ints$dist != 0,]
    }
	
    if(trans==TRUE)
    {
        for(chr in chrs)
        {
            print(chr)
   	    d1 <- set1[as.character(set1$chrom) == chr,]
	    d2 <- set2[as.character(set2$chrom) != chr,]
		
	    if(nrow(d1) == 0 | nrow(d2) == 0){next}
			
	    d_inds = as.data.frame(expand.grid(1:nrow(d1), 1:nrow(d2)))
	    d_ints_for = data.frame(chrom1=d1[d_inds[,1],'chrom'], start1=d1[d_inds[,1],'start'], end1=d1[d_inds[,1],'end'], chrom2=d2[d_inds[,2],'chrom'], start2=d2[d_inds[,2],'start'], end2=d2[d_inds[,2],'end'])
	    d_ints_rev = data.frame(chrom1=d2[d_inds[,2],'chrom'], start1=d2[d_inds[,2],'start'], end1=d2[d_inds[,2],'end'], chrom2=d1[d_inds[,1],'chrom'], start2=d1[d_inds[,1],'start'], end2=d1[d_inds[,1],'end'])
	    d_ints <- rbind(d_ints_for,d_ints_rev)
	    d_ints <- d_ints[!duplicated(paste(d_ints$chrom1,d_ints$start1,d_ints$end1,d_ints$chrom2,d_ints$start2,d_ints$end2,sep='_')),]
	    d_ints$dist <- NA
			
	    print(head(d_ints))
			
	    if(is.null(ints))
	    {
	        ints <- d_ints
	    }else{
	        ints <- rbind(ints,d_ints)
	        ints <- ints[!duplicated(paste(ints$chrom1,ints$start1,ints$end1,ints$chrom2,ints$start2,ints$end2,sep='_')),]
	    }
	}
    }
	
    if(retNames == TRUE)
    {
		
        if(names1[1]=='')
        {
            names1 <- paste(paste0('p',1:nrow(set1)),set1$chrom,set1$start,set1$end,sep='_')
        }
        if(names2[1]=='')
        {
            names2 <- paste(paste0('p',1:nrow(set2)),set2$chrom,set2$start,set2$end,sep='_')
        }
		
        outNames1 <- c()
        outNames2 <- c()
        for (i in 1:nrow(ints))
        {
    	    x <- ints[i,]
	    s1 <- paste(x['start1'],x['end1'],sep='_')
	    s2 <- paste(x['start2'],x['end2'],sep='_')
	    oN1 <- unlist(strsplit(names1[grep(s1,names1)],'_'))[1]
	    oN2 <- unlist(strsplit(names2[grep(s2,names2)],'_'))[1]
	    outNames1 <- append(outNames1,oN1)
	    outNames2 <- append(outNames2,oN2)
	}
	
	ints$name1 <- outNames1
	ints$name2 <- outNames2
    }

    return(ints)
}
### END : build2D ###

### BEGIN : plotRectangles2D ###
plotRectangles2D <- function(data, int2d=limits2d, fc=0.5, col='yellow', lwd=0.5, lty=1, image=TRUE, side='both', shiftEnd=0, shiftStart=0,shade=FALSE)
{
	
    if(is.null(data)){return()}
	
    #### Normalize @ 1 by 1 for image...
    xLims <- c(int2d[2:3])
    xL <- as.numeric(as.character(int2d[3] - int2d[2]))
    yLims <- c(int2d[5:6])
    yL <- as.numeric(as.character(int2d[6] - int2d[5]))
	
    start <- as.numeric(as.character(xLims[1]))
    end <- as.numeric(as.character(xLims[2]))
	
    for(g in 1:nrow(data))
    {
 	row <- data[g,]

 	if(image==FALSE)
	{
  	    x1 <- row$start1
	    y1 <- row$end1
	    x2 <- row$start2
	    y2 <- row$end2
	}else if (image==TRUE){
	    x1 <- (row$start1 - as.numeric(as.character(xLims[1])))/xL
	    y1 <- (row$end1 - as.numeric(as.character(xLims[1])))/xL
	    x2 <- (row$start2 - as.numeric(as.character(yLims[1])))/yL
	    y2 <- (row$end2 - as.numeric(as.character(yLims[1])))/yL
	}else if (image=='TRI'){
			
	    if(row$start1 != row$start2)
	    {
	    	x1 <- row$start1 + (row$start2-row$start1)/2
		x2 <- row$start1 + (row$end2-row$start1)/2
		x3 <- row$end1 + (row$end2-row$end1)/2
		x4 <- row$end1 + (row$start2-row$end1)/2
				
		y1 <- row$start2-row$start1
		y2 <- row$end2-row$start1
		y3 <- row$end2-row$end1
		y4 <- row$start2-row$end1
				
		segments(x1,y1,x2,y2,col=col,lwd=lwd,lty=lty)
		segments(x2,y2,x3,y3,col=col,lwd=lwd,lty=lty)
			
		segments(x3,y3,x4,y4,col=col,lwd=lwd,lty=lty)
		segments(x4,y4,x1,y1,col=col,lwd=lwd,lty=lty)
				
		if(shade == TRUE)
		{
		    print('filled...')
		    polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=1)
		}
	    }else{
		x1 <- row$start1
		x2 <- row$end2
		mid <- x1 + (x2-x1)/2
		y1 <- 0
		y2 <- (x2-x1)
				
		if(nrow(data) > 1)
		{
		    if(g == 1 & row$start1 < start)
		    {
			x1 <- row$end2
			x2 <- row$end2
			mid <- start + (x2-start)/2
			y1 <- 0
			y2 <- (x2-start)
		    }else if(g == nrow(data) & row$end2 > end){
			x1 <- row$start2
			x2 <- row$start2
			mid <- x1 + (end-x1)/2
			y1 <- 0
			y2 <- (end-x1)
  		    }
		}
				
		if(side=='top')
		{
		    segments(x1,y1,mid,y2,col=col,lwd=lwd,lty=lty)
		    segments(mid,y2,x2,y1,col=col,lwd=lwd,lty=lty)
		}else if(side=='bottom'){
		    segments(x2,y1,mid,-y2,col=col,lwd=lwd,lty=lty)
		    segments(mid,-y2,x1,y1,col=col,lwd=lwd,lty=lty)
		}else{
		    segments(x1,y1,mid,y2,col=col,lwd=lwd,lty=lty)
		    segments(mid,y2,x2,y1,col=col,lwd=lwd,lty=lty)
		    segments(x2,y1,mid,-y2,col=col,lwd=lwd,lty=lty)
		    segments(mid,-y2,x1,y1,col=col,lwd=lwd,lty=lty)
		}
	    }
	}
		
	if (image!='TRI')
	{
	    if(side=='top')
	    {
		segments(x1,x2,x1,y2,col=col,lwd=lwd,lty=lty)
		segments(x1,y2,y1,y2,col=col,lwd=lwd,lty=lty)
	    }else if(side=='bottom'){
	     	segments(y1,y2,y1,x2,col=col,lwd=lwd,lty=lty)
		segments(x1,x2,y1,x2,col=col,lwd=lwd,lty=lty)
	    }else{
		segments(x1,x2,x1,y2,col=col,lwd=lwd,lty=lty)
		segments(x1,y2,y1,y2,col=col,lwd=lwd,lty=lty)
		segments(y1,y2,y1,x2,col=col,lwd=lwd,lty=lty)
		segments(x1,x2,y1,x2,col=col,lwd=lwd,lty=lty)
	    }
	}
    }
}
### END : plotRectangles2D ###

### BEGIN : expand1D ###
expand1D <- function(set,expand,center=1)
{
    require(plyr)
    set$dist <- set$end - set$start
	
    set$mid <- set$start + round_any((set$end-set$start)/2,center)
	
    set <- set[order(set$mid),]
	
    mids <- set$mid
    min <- 10e10
	
    if(nrow(set) > 1)
    {
	for (i in 1:(length(mids)-1))
	{
	    dist <- mids[i+1] - mids[i]
	    if(dist < min ){min <- dist}
	}
		
    	if(min < expand*2)
	{
 	    print(paste('************************************************',min))
        }
    }
	
    set$start <- set$mid - expand
    set$end <- set$mid + expand
	
    set <- set[set$start > 0,]
	
    for (chr in as.character(gintervals.all()$chrom))
    {
	maxSize <- as.numeric(as.character(gintervals.all()[gintervals.all()$chrom == chr,'end']))
	if(length( which(set$chrom == chr & set$end > maxSize) ) > 0)
	{
	    set <- set[-which(set$chrom == chr & set$end > maxSize),]
	}
    }
	
    set <- gintervals.canonic(set)
	
    return(set)
}
### END : expand1D ###

### BEGIN : expand2D ###
expand2D <- function(set,expand,selectOn=NULL,higher=TRUE)
{

    if(nrow(set) == 0)
    {
  	return(NULL)
    }

    set[is.na(set)] <- 0
	
    set$start1 <- set$start1 - expand
    set$end1 <- set$end1 + expand
    set$start2 <- set$start2 - expand
    set$end2 <- set$end2 + expand
	
    if(set$start1 < 0){set$start1 <- 0}
    if(set$start2 < 0){set$start2 <- 0}
	
    set <- set[set$start1 >= 0,]
	
    for (chr in as.character(gintervals.all()$chrom))
    {
    	maxSize <- as.numeric(as.character(gintervals.all()[gintervals.all()$chrom == chr,'end']))
		
	set[which(set$chrom1 == chr & set$end1 > maxSize),'end1'] <- maxSize
	set[which(set$chrom1 == chr & set$end2 > maxSize),'end2'] <- maxSize
		
		
    }
	
    set$mid1 <- set$start1 + round_any((set$end1-set$start1)/2,1)
    set$mid2 <- set$start2 + round_any((set$end2-set$start2)/2,1)
	
    if(nrow(set) == 1)
    {
	return(set)
    }
	
    print('starting')
    print(nrow(set))
	
    done <- 1
    while(done > 0)
    {
    	set <- set[order(set$chrom1,set$mid1+set$mid2),]
 	print(nrow(set))
	remove <- c()
	d1 <- 0
	d2 <- 0
		
	if(nrow(set) < 2){done <- 0; break}
		
	for (i in 1:(nrow(set)-1))
	{
	    if(set[i,'chrom1'] != set[i+1,'chrom1']){next}
			
  	    d1 <- set[i+1,'mid1'] - set[i,'mid1']
	    d2 <- set[i+1,'mid2'] - set[i,'mid2']
			
	    if(d1 <= expand*2 & d2 <= expand*2)
	    {
	    	if(!is.null(selectOn))
		{
		    if(set[i,which(colnames(set)==selectOn)] >= set[i+1,which(colnames(set)==selectOn)])
		    {
		    	if(higher == TRUE)
			{
			    remove <- append(remove,i+1)
			}else{
			    remove <- append(remove,i)
			}
		    }else{
		        if(higher == TRUE)
		        {
		            remove <- append(remove,i)
		        }else{
			    remove <- append(remove,i+1)
			}
		    }
		}else{
		    remove <- append(remove,i)
		}
	    }
	}
		
	remove <- unique(remove)
		
	set$d1 <- d1
	set$d2 <- d2
	done <- length(remove)
		
	if(done != 0)
	{
	    set <- set[-remove,]
	}
    }
	
    done <- 1
    while(done > 0)
    {
		
        set <- set[order(set$chrom1,set$mid2,set$mid1),]
	print(nrow(set))
	remove <- c()
	d1 <- 0
	d2 <- 0
		
	if(nrow(set) < 2){done <- 0; break}
		
	for (i in 1:(nrow(set)-1))
	{
	    if(set[i,'chrom1'] != set[i+1,'chrom1']){next}
		
  	    d1 <- set[i+1,'mid1'] - set[i,'mid1']
	    d2 <- set[i+1,'mid2'] - set[i,'mid2']
			
	    if(d1 <= expand*2 & d2 <= expand*2)
	    {
	    	if(!is.null(selectOn))
		{
		    if(set[i,which(colnames(set)==selectOn)] >= set[i+1,which(colnames(set)==selectOn)])
		    {
		    	if(higher == TRUE)
			{
			    remove <- append(remove,i+1)
			}else{
			    remove <- append(remove,i)
			}
		    }else{
			if(higher == TRUE)
		        {
			    remove <- append(remove,i)
			}else{
			    remove <- append(remove,i+1)}
			}
		}else{
		    remove <- append(remove,i)
 	        }
	    }
        }
	remove <- unique(remove)
		
	set$d1 <- d1
	set$d2 <- d2
	done <- length(remove)
		
	if(done != 0)
	{
	    set <- set[-remove,]
	}
    }
	
    return(set)
}
### END : expand2D ###

### Insulation ###
### BEGIN : gtrack.2d.gen_insu_track ###
gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description="")
{
    k_reg = 10
    if (length(track_nm) > 1)
    {
 	prof = gtrack.2d.gen_joint_insu_prof(track_nm, scale, res, min_diag_d)
    }else{
	print("single source track")
	prof = gtrack.2d.gen_insu_prof(track_nm, scale, res, min_diag_d)
    }
	
    message("names ", paste(names(prof), collapse=","))
    names(prof)[1] = "chrom"
    names(prof)[2] = "start"
    names(prof)[3] = "end"
	
    gtrack.create_sparse(track=new_track, prof[,c(1,2,3)], log2(prof$obs_ins/(prof$obs_big+k_reg)), description=description)
}
### END : gtrack.2d.gen_insu_track ###

### BEGIN : gtrack.2d.gen_joint_insu_prof ###
gtrack.2d.gen_joint_insu_prof = function(track_nms, scale, res, min_diag_d = 1000)
{
    ins = gtrack.2d.gen_insu_prof(track_nms[1], scale, res, min_diag_d)
	
    for (track in track_nms[-1])
    {
	t_ins =  gtrack.2d.gen_insu_prof(track, scale, res, min_diag_d)
	ins$obs_big = rowSums(cbind(ins$obs_big, t_ins$obs_big), na.rm=T)
	ins$obs_ins = rowSums(cbind(ins$obs_ins, t_ins$obs_ins), na.rm=T)
    }
	
    ins$obs_big[ins$obs_big == 0] = NA
    ins$obs_ins[ins$obs_ins == 0] = NA
    return(ins)
}
### END : gtrack.2d.gen_joint_insu_prof ###

### BEGIN : gtrack.2d.gen_insu_prof ###
gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000)
{
    iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)
	
    iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end,
  			    chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
    if(length(gvtrack.ls("obs_big")) == 1){gvtrack.rm("obs_big")}
	
    if(length(gvtrack.ls("obs_ins")) == 1){gvtrack.rm("obs_ins")}
	
    gvtrack.create("obs_big", track_nm, "weighted.sum")
    gvtrack.create("obs_ins", track_nm, "weighted.sum")
    gvtrack.iterator.2d("obs_big", 
	eshift1=scale, sshift1=-scale, 
	eshift2=scale, sshift2=-scale)
    gvtrack.iterator.2d("obs_ins", 
	eshift1=0, sshift1=-scale, 
	eshift2=scale, sshift2=0)
	
    message("will iter on ", dim(iter_2d)[1])
    ins = gextract("obs_big", "obs_ins", gintervals.2d.all(), iterator=iter_2d, band=c(-scale*2,0))
    ins_diag = gextract("obs_big", "obs_ins", gintervals.2d.all(), iterator=iter_2d, band=c(-min_diag_d,0))
	
    ###fixing NA bug
    ins[is.na(ins)] = 0
    ins_diag[is.na(ins_diag)] = 0
    ###
    ins$obs_big = ins$obs_big - ins_diag$obs_big
    ins$obs_ins = ins$obs_ins - ins_diag$obs_ins
    message("will retrun ins with ", dim(ins)[1], " rows")
    return(ins)
}
### END : gtrack.2d.gen_insu_prof ###
