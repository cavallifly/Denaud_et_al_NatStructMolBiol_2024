require(misha)
require(shaman)

require(doParallel)

##########################################
### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('./scripts/auxFunctions.R')
options(scipen=20,gmax.data.size=0.5e8)

refData <- list(
	hic_larvae_DWT_merge_dm6_BeS   = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1'),
	hic_LD_gypsy1_merge_dm6_BeS    = c('hic.hic_LD_gypsy1_Rep1_dm6_BeS.hic_LD_gypsy1_Rep1_dm6_BeS_1','hic.hic_LD_gypsy1_Rep2_dm6_BeS.hic_LD_gypsy1_Rep2_dm6_BeS_1','hic.hic_LD_gypsy1_Rep3_dm6_BeS.hic_LD_gypsy1_Rep3_dm6_BeS_1'),
        hic_LD_gypsy2Enh_merge_dm6_BeS = c('hic.hic_LD_gypsy2Enh_Rep1_dm6_BeS.hic_LD_gypsy2Enh_Rep1_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep3_dm6_BeS.hic_LD_gypsy2Enh_Rep3_dm6_BeS_1'),
        hic_LD_gypsy3_merge_dm6_SD     = c('hic.hic_LD_gypsy3_Rep1_dm6_SD.hic_LD_gypsy3_Rep1_dm6_SD_1','hic.hic_LD_gypsy3_Rep2_dm6_SD.hic_LD_gypsy3_Rep2_dm6_SD_1')
	)

obsData <- list(
	hic_larvae_DWT_merge_dm6_BeS   = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1')
	)

print(refData)
print(obsData)
refTracks  <- refData
obsTracks <- obsData
print(refTracks)
print(obsTracks)

registerDoParallel(cores=12)

chr <- 'chr2L'
targetRegion2D <- g2d(gintervals(chr,0,gintervals.all()[gintervals.all()$chrom == chr,'end']))
print(targetRegion2D)

sampleCoverageAll <- list()
sampleInts        <- list()
targetInts <- list()
refTrackNames <- list()

for(refTrack in names(refTracks))
{
    sampleCoverage <- 0
    print(refTrack)
    print(refTracks[[refTrack]])
    print(typeof(refTracks[[refTrack]]))
    refTrackName  <- gsub('_merge_mm10_IJ','',gsub('chic_','',refTrack))
    refTrackNames[[refTrack]] <- refTrackName

    for(currentTrack in refTracks[[refTrack]])
    {
        splitName <- unlist(strsplit(currentTrack,'\\.'))
	setName <- splitName[length(splitName)]
        print(setName)
        setName <- gsub('_merge_dm6_BS','',gsub('chic_','',setName))
        print(setName)

	targetData <- gextract(currentTrack,targetRegion2D,colnames=("obs"))
        targetInts[[setName]] <- sum(targetData$obs)/2
        print(sum(targetData$obs)/2)
        sampleCoverage <- sampleCoverage + sum(targetData$obs)/2
    }

    sampleInts[[refTrack]] <- sampleCoverage
    print(sampleInts[[refTrack]])
}
print(sampleInts)

minInts <- min(unlist(sampleInts))
minSample <- names(sampleInts)[which(unlist(sampleInts) == min(unlist(sampleInts)))]
print(minSample)

sampleRatios <- list()
for(set in names(refTracks))
{
    if(set == minSample)
    {
	sampleRatio <- 1
    }else{
	sampleRatio <- sampleInts[[minSample]]/sampleInts[[set]]
    }
 	
    print(set)
    print(sampleRatio)

    sampleRatios[[set]] <- sampleRatio

    outTrack <- paste0('hic.diffMaps.',set,'_randMin')
    print(outTrack)

    if(gtrack.exists(outTrack)){next}
 	
    for(currentTrack in refTracks[[set]])
    {
 	
 	setName <- unlist(strsplit(currentTrack,'\\.'))[3]
 	print(paste('-----',setName))
	
 	outTXT <- paste0('scoreTemp/',set,'_rep',which(refTracks[[set]] == currentTrack),'_minData.txt')		
 	if(file.exists(outTXT)){print('done'); next}
 		
 	targetData <- gextract(currentTrack,targetRegion2D)
 	targetData <- gintervals.canonic(targetData)
	print(head(targetData))

	targetData <- targetData[targetData$start2 > targetData$start1,]
 	targetData$obs <- 1
	print(nrow(targetData))
	print(nrow(targetData)*sampleRatio)
	minData <- targetData[sample(nrow(targetData), nrow(targetData)*sampleRatio),1:7]
	colnames(minData)[7] <- 'obs'
 	options(scipen=1)
 	write.table(minData,file=outTXT,quote=F, row.names=F, sep="\t")
    }
 	
    dataFiles <- list.files('scoreTemp/',pattern=paste0(set,'_rep'),full.names=TRUE)
 	
    if(gtrack.exists(outTrack)){gtrack.rm(outTrack, force=TRUE)}
    gtrack.2d.import_contacts(outTrack, paste("Random minimum 2d track for", set), dataFiles)
}
 
randInts <- list()
for(currentTrack in gtrack.ls('hic.diffMaps'))
{
    if(length(grep('randMin',currentTrack)) == 0){next}
    setName <- unlist(strsplit(currentTrack,'\\.'))[3]
    print(setName)
    gvtrack.create('vTrack',currentTrack,'area')
    targetData <- gextract('vTrack',targetRegion2D,iterator=targetRegion2D)
    randInts[[setName]] <- targetData$vTrack/2
}
print(randInts)

pdf(file='sandrine_counts.pdf',width=14,height=11)
par(mfrow=c(2,2),mai=rep(c(1.65,0.75),2))
barplot(unlist(targetInts),las=2,main=paste0("target region (",chr,") ints"), col=rep(c('darkslategrey','violet','lavender'),length(names(refTrackNames))))

yMax <- max(unlist(sampleInts))*1.2
barplot(unlist(sampleInts),beside=T,las=2,main=paste0("target region (",chr,") ints"), ylim=c(0,yMax),space=0,cex.main=1, cex.names=1.5,las=2,col=colorRampPalette(c('darkslategrey','violet','lavender'))(length(names(refTrackNames))))
cnts <- paste0(signif(unlist(sampleInts),3)/1e6,'')
sapply(1:length(names(sampleInts)),function(y){text(x=y-0.5,y=sampleInts[[y]],labels=cnts[y],pos=3,cex=1.5)})

yMax <- max(unlist(sampleRatios))*1.2
barplot(unlist(sampleRatios),beside=T,las=2,main=paste0('target region (',chr,') ints'), ylim=c(0,yMax),space=0,cex.main=1, cex.names=1.5,las=2,col=colorRampPalette(c('darkslategrey','violet','lavender'))(length(names(sampleRatios))))
cnts <- paste0(signif(unlist(sampleRatios),3)*100,'%')
sapply(1:length(names(sampleRatios)),function(y){text(x=y-0.5,y=sampleRatios[[y]],labels=cnts[y],pos=3,cex=1.5)})
 
yMax <- max(unlist(randInts))*1.2
barplot(unlist(randInts),beside=T,las=2,main=paste0('target region (',chr,') DOWN-SAMPLED ints'), ylim=c(0,yMax),space=0,cex.main=1, cex.names=1.5,las=2,col=colorRampPalette(c('darkslategrey','violet','lavender'))(length(names(randInts))))
cnts <- paste0(signif(unlist(randInts),3)/1e6,'')
sapply(1:length(names(randInts)),function(y){text(x=y-0.5,y=randInts[[y]],labels=cnts[y],pos=3,cex=1.5)})
dev.off()

step <- 1e6
expand <- step/3
k    <- 250
kexp <- NA
 
registerDoParallel(cores=12)

randMinTracks <- gtrack.ls('randMin')
print(randMinTracks)

chr   <- 'chr2L'
start <- 10e6
end   <- 20e6

chrIntrv <- gintervals(chr,start,end)
chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step))
print(chrIntrv)

#print("Computing differential maps!")

for(obsTrack in randMinTracks)
{
    print(obsTrack)
    for(expTrack in randMinTracks)
    {
        print(expTrack)

	if(obsTrack == expTrack){next}	

	hicRoot <- "hic.diffMaps."
 		
	set1 <- gsub('_randMin','',unlist(strsplit(obsTrack,'\\.'))[3])
	set2 <- gsub('_randMin','',unlist(strsplit(expTrack,'\\.'))[3])

	trackName <- paste0(set1,'_vs_',set2,'_',chr,'_score_k',k,'_kexp',kexp,"_",step/1e6,'Mb')
	print(trackName)
	outTrack <- paste0(hicRoot,trackName)		
	print(outTrack)
	hicName <- paste0(set1,'_vs_',set2)
	work_dir <- paste0("./_tmp_",trackName,"/")

	if(!dir.exists(work_dir)){dir.create(work_dir, mode="7777", recursive=TRUE);}

	if(gtrack.exists(outTrack))
	{
	    print(paste0(outTrack," exists"))
	    next
	}

	foreach(i=1:nrow(chrIter)) %dopar% {
	int1 <- chrIter[i,]
	int2 <- expand2D(int1,expand)[,1:6]
			
	outFile <- paste0(work_dir,chr,'_',i,'.scores')
	if(file.exists(outFile))
	{
	    print('done...')
   	    return(i)
	}
 			
	print('scoring...')
	scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2, k=k, k_exp=kexp)
	data <- scores$points
	write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)
 			
	return(i)
    }
	
    files <- list.files(work_dir, full.names=T,pattern='scores')
 		
    ### Single Import
    gtrack.2d.import(outTrack, paste("Score track for", hicName,' - ',chr), files)
    ### Add Int AND Reverse INT
    # gtrack.2d.import_contacts(outTrack, paste("import diffMaps scores for ", trackName), files)

    ### Remove temp files....
    for (f in files){file.remove(f)}
    }
}
quit()
