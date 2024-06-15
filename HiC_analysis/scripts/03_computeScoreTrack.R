require(misha)
require(shaman)

require(doParallel)

args = commandArgs(trailingOnly=TRUE)

##########################################
### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('./scripts/auxFunctions.R')
options(scipen=20,gmax.data.size=0.6e8,shaman.sge_support=1)

#chr <- 'chr2L'
#start <- 10e6
#end <- 20e6

step   <- 1000000 # Size of the chunk to divide the calculation
expand <- step/3  # Size of skin around the main chunk to look for expected counts v2
k      <- 250     # Number of neighbours to look for
k_exp  <- 2*k 

samples = c("larvae_DWT","LD_DRE","LD_gypsy1","LD_gypsy2Enh","pupae_DWT","LE_WT","LD_gypsy3")
print(samples)	

for(sample in samples)
{
    print(sample)	
    if(sample == "larvae_DWT")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep1_dm6_BeS/hic_larvae_DWT_Rep1_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep2_dm6_BeS/hic_larvae_DWT_Rep2_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep3_dm6_BeS/hic_larvae_DWT_Rep3_dm6_BeS_1.track:
	nreplicates = 3
	author      = "BeS"
    }

    if(sample == "LD_DRE")
    {
    	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_DRE_Rep1_dm6_BeS/hic_LD_DRE_Rep1_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_DRE_Rep2_dm6_BeS/hic_LD_DRE_Rep2_dm6_BeS_1.track:
	nreplicates = 2
	author      = "BeS"
    }

    if(sample == "LD_gypsy1")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy1_Rep1_dm6_BeS/hic_LD_gypsy1_Rep1_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy1_Rep2_dm6_BeS/hic_LD_gypsy1_Rep2_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy1_Rep3_dm6_BeS/hic_LD_gypsy1_Rep3_dm6_BeS_1.track:
	nreplicates  = 3
	author       = "BeS"
    }

    if(sample == "LD_gypsy2Enh")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy2Enh_Rep1_dm6_BeS/hic_LD_gypsy2Enh_Rep1_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy2Enh_Rep2_dm6_BeS/hic_LD_gypsy2Enh_Rep2_dm6_BeS_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy2Enh_Rep3_dm6_BeS/hic_LD_gypsy2Enh_Rep3_dm6_BeS_1.track:
	nreplicates  = 3
	author       = "BeS"
    }

    if(sample == "pupae_DWT")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_pupae_DWT_Rep1_dm6_BeS/hic_pupae_DWT_Rep1_dm6_BeS_1.track
	#./mishaDB/trackdb/dm6/tracks/hic/hic_pupae_DWT_Rep2_dm6_BeS/hic_pupae_DWT_Rep2_dm6_BeS_1.track
	#./mishaDB/trackdb/dm6/tracks/hic/hic_pupae_DWT_Rep3_dm6_BeS/hic_pupae_DWT_Rep2_dm6_BeS_1.track
	nreplicates  = 3
	author       = "BeS"
    }

    if(sample == "LE_WT")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LE_WT_Rep1_dm6_YO/hic_LE_WT_Rep1_dm6_YO_1.track
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LE_WT_Rep2_dm6_YO/hic_LE_WT_Rep2_dm6_YO_1.track
    	#./mishaDB/trackdb/dm6/tracks/hic/hic_LE_WT_Rep3_dm6_VL/hic_LE_WT_Rep3_dm6_VL_1.track
    	nreplicates  = 3
    	author       = "YOVL"
    }

    if(sample == "LD_gypsy3")
    {
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy3_Rep1_dm6_SD/hic_LD_gypsy3_Rep1_dm6_SD_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy3_Rep1res_dm6_SD/hic_LD_gypsy3_Rep2_dm6_SD_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy3_Rep2_dm6_SD/hic_LD_gypsy3_Rep1_dm6_SD_1.track:
	#./mishaDB/trackdb/dm6/tracks/hic/hic_LD_gypsy3_Rep2res_dm6_SD/hic_LD_gypsy3_Rep2_dm6_SD_1.track:
	nreplicates  = 4
	author       = "SD"
    }

    for(chr in gintervals.all()$chrom)
    #for(chr in c(chr))
    {
        print(chr)

        chrIntrv <- gintervals.all()[gintervals.all()$chrom == chr,]
        chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step)) #,band=c(-10e10,0))
        print(chr);print(chrIntrv$start);print(chrIntrv$end)
        print(nrow(chrIter))

	### SCORE MAPS
	if(nreplicates == 2)
    	{
	    obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"))
	}
	if(nreplicates == 3)
	{        
            obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep3_dm6_",author,".hic_",sample,"_Rep3_dm6_",author,"_1"))
	    if(sample == "LE_WT")
	    {
	        obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_YO.hic_",sample,"_Rep1_dm6_YO_1"),paste0("hic.hic_",sample,"_Rep2_dm6_YO.hic_",sample,"_Rep2_dm6_YO_1"),paste0("hic.hic_",sample,"_Rep3_dm6_VL.hic_",sample,"_Rep3_dm6_VL_1"))
	    }
        }
        if(nreplicates == 4)
        {
            obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep1res_dm6_",author,".hic_",sample,"_Rep1res_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2res_dm6_",author,".hic_",sample,"_Rep2res_dm6_",author,"_1"))
        }

        setName  <- paste0("hicScores_",sample,"_merge_dm6_",author)
        if(sample == "LE_WT")
	{
	    if(author == "YOVL")
	    {
                setName  <- paste0("hicScores_",sample,"_merge_dm6_YOVL")
	    }
        }

	expTrack <- paste0(obsTrack,"_shuffle_500Small_1000")
	print(obsTrack)
	print(setName)
	print(expTrack)

	hicName <- paste0("hic")
	
	trackName <- paste0(setName,"_",chr,"_",chrIntrv$start/1e6,"_",chrIntrv$end/1e6,"Mb_score_k",k,"_kexp",k_exp,"_step",step/1e6,"Mb")
	outTrack <- paste0(hicName,".",setName,".",trackName)
	print(outTrack)

	work_dir <- paste0("./_tmp_",trackName,"/")
	if(!dir.exists(work_dir)){dir.create(work_dir, mode="7777", recursive=TRUE);}

	if(gtrack.exists(outTrack))
	{
            print(paste0(outTrack," exists!"))
            next
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
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k,k_exp=k_exp)
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
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k,k_exp=k_exp)
	    data <- scores$points
	    write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	    return(i)
	}

	files <- list.files(work_dir, full.names=T,pattern="scores")
	
	trackFolder <- paste0(dbDir,"tracks/",gsub("\\.","/",hicName),"/",setName,"/")
	if(!dir.exists(trackFolder)){dir.create(trackFolder, mode="7777", recursive=TRUE);}

	### Single Import
	gtrack.2d.import(outTrack, paste("Score track for", hicName," - ",chr), files)
	### Add Int AND Reverse INT
	# gtrack.2d.import_contacts(outTrack, paste("import diffMaps scores for ", trackName), files)

	### Remove temp files....
	for (f in files){file.remove(f)}
    }
}
