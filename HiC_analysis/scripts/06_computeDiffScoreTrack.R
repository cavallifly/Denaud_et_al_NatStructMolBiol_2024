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
        hic_LD_gypsy1_merge_dm6_BeS    = c('hic.hic_LD_gypsy1_Rep1_dm6_BeS.hic_LD_gypsy1_Rep1_dm6_BeS_1','hic.hic_LD_gypsy1_Rep2_dm6_BeS.hic_LD_gypsy1_Rep2_dm6_BeS_1','hic.hic_LD_gypsy1_Rep3_dm6_BeS.hic_LD_gypsy1_Rep3_dm6_BeS_1'),
	hic_LD_gypsy2Enh_merge_dm6_BeS = c('hic.hic_LD_gypsy2Enh_Rep1_dm6_BeS.hic_LD_gypsy2Enh_Rep1_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep3_dm6_BeS.hic_LD_gypsy2Enh_Rep3_dm6_BeS_1'),
        hic_LD_gypsy3_merge_dm6_SD     = c('hic.hic_LD_gypsy3_Rep1_dm6_SD.hic_LD_gypsy3_Rep1_dm6_SD_1','hic.hic_LD_gypsy3_Rep2_dm6_SD.hic_LD_gypsy3_Rep2_dm6_SD_1')
	)

obsData <- list(
	hic_larvae_DWT_merge_dm6_BeS   = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1')
	)

print(refData)
print(obsData)
refTracks <- refData
obsTracks <- obsData
print(refTracks)
print(obsTracks)

### VS MAPS variable K
registerDoParallel(cores=12)
step   <- 3e6
expand <- step/3
k      <- 250
kexp  <- k # kexp=k has to be used to compare datasets

geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
rownames(geneCoordinates) <- geneCoordinates$geneName
print(geneCoordinates['dac',1:3])

chr <- 'chr2L'
start <- 16354000
end   <- 16500000

chrIntrv <- gintervals(chr,start,end)
chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step)) #,band=c(-10e10,0))
print(dim(chrIter))

sampleCoverageAll <- list()
for(refTrack in names(refTracks))
{
    print(refTrack)
    print(refTracks[[refTrack]])    
    print(typeof(refTracks[[refTrack]]))
    targetRegion2D <- g2d(chrIntrv)
    print(targetRegion2D)

    tracks = refTracks[[refTrack]]
    vTracks <- paste0('v_',tracks,'_all')
    print(vTracks)
    for(j in 1:length(tracks))
    {
       #vTrack <- gvtrack.create(vTracks[j],tracks[j],"area")
       vTrack <- gvtrack.create(paste0("ref",j),tracks[j],"area")       
       vTracks[j] <- paste0("ref",j)
    }
    targetData <- gextract(vTracks,targetRegion2D,iterator=targetRegion2D)
    targetData <- targetData[vTracks]
    print(head(targetData))
    sampleCoverageAll[[refTrack]] <- 2*rowSums(targetData)
}
print(sampleCoverageAll)

allTracks <- gtrack.ls('hic.diffMaps')
for(obsTrack in names(obsTracks))
{
    for(expTrack in names(refTracks))
    {	
		
        hicName <- 'hic'
	print(obsTrack)
	print(expTrack)	
	set1 <- obsTrack
	set2 <- expTrack	
	#set1 <- unlist(strsplit(obsTrack,'\\.'))[2]
        #set2 <- unlist(strsplit(expTrack,'\\.'))[2]
	if(set1 == set2){next}
	trackName <- paste0(set1,'_vs_',set2,'_',chr,'_score_k',k,'_kexp',kexp,'_3Mb')
        outTrack <- paste0(hicName,'.diffMaps.',trackName)
	print(outTrack)

        if(gtrack.exists(outTrack)){next}

        work_dir <- paste0('_tmp_',outTrack,'/')
	print(obsTracks[[obsTrack]])
	print(refTracks[[expTrack]])

	if(!dir.exists(work_dir)){dir.create(work_dir, mode="7777", recursive=TRUE);}

	if(sampleCoverageAll[[obsTrack]] < sampleCoverageAll[[expTrack]])
	{
	    print(paste0("RANDOM SAMPLING ",obsTrack)
            sampleRatio <- sampleCoverageAll[[obsTrack]]/sampleCoverageAll[[expTrack]]
            print(sampleRatio)
            obsTrack <- 
            expTrack <- expTrack
        }else{
            print("RANDOM SAMPLING ",expTrack)
            sampleRatio <- sampleCoverageAll[[expTrack]]/sampleCoverageAll[[obsTrack]]
            print(sampleRatio)
            expTrack <- obsTracks[grep(paste0(setName2,'_rep'),obsTracks)]
            sampleTracks <- wtTracks[grep(stage,wtTracks)]
        }

        registerDoParallel(cores=5)

        foreach(i=1:nrow(chrIter)) %dopar% {
            int1 <- chrIter[i,]	
	    int2 <- expand2D(int1,expand)[,1:6]
	    print(paste0("Interval to score"))
	    print(int1)	
	    print(paste0("Interval + skin"))	
	    print(int2)

	    outFile <- paste0(work_dir,chr,"_",i,".scores")
	    if(file.exists(outFile))
	    {
	        print("done...")
		return(i)
	    }

  	    genDist <- abs(int1$start2-int1$start1)
	    if(genDist <= 500e3)
	    {
	        print(paste0("Closer than ",500e3,"bp. Leave it for later!"))	
		return(i)
	    }

	    print(paste0("scoring portion ",i,"..."))
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTracks[[obsTrack]],exp_track_nms=refTracks[[expTrack]],focus_interval=int1,regional_interval=int2,k=k,k_exp=kexp)
	    data <- scores$points
	    write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	    return(i)
	}

        registerDoParallel(cores=1)	

        foreach(i=1:nrow(chrIter)) %dopar% {
            int1 <- chrIter[i,]	
	    int2 <- expand2D(int1,expand)[,1:6]
	    print(paste0("Interval to score"))
	    print(int1)	
	    print(paste0("Interval + skin"))	
	    print(int2)

	    outFile <- paste0(work_dir,chr,"_",i,".scores")
	    if(file.exists(outFile))
	    {
	        print("done...")
		return(i)
	    }

	    genDist <- abs(int1$start2-int1$start1)
	    if(genDist > 500e3)
	    {
	        print(paste0("Further than ",500e3,"bp. Already done!"))
		return(i)
	    }

	    print(paste0("scoring portion ",i,"..."))
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTracks[[obsTrack]],exp_track_nms=refTracks[[expTrack]],focus_interval=int1,regional_interval=int2,k=k,k_exp=kexp)
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
