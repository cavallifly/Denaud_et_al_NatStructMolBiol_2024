require(misha)
require(shaman)

require(doParallel)

########################
### misha working DB ###
mDBloc <-  './mishaDB/trackdb/'
db     <- 'dm6'
dbDir  <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('./scripts/auxFunctions.R')
options(scipen=20,gmax.data.size=1.e9)

refData <- list(
	hic_LE_WT        = c("hic.hic_LE_WT_Rep1_dm6_YO.hic_LE_WT_Rep1_dm6_YO_1","hic.hic_LE_WT_Rep2_dm6_YO.hic_LE_WT_Rep2_dm6_YO_1","hic.hic_LE_WT_Rep3_dm6_VL.hic_LE_WT_Rep3_dm6_VL_1"),
	hic_larvae_DWT   = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1'),
	hic_pupae_DWT    = c('hic.hic_pupae_DWT_Rep1_dm6_BeS.hic_pupae_DWT_Rep1_dm6_BeS_1','hic.hic_pupae_DWT_Rep2_dm6_BeS.hic_pupae_DWT_Rep2_dm6_BeS_1','hic.hic_pupae_DWT_Rep3_dm6_BeS.hic_pupae_DWT_Rep3_dm6_BeS_1'),
	hic_LD_DRE       = c("hic.hic_LD_DRE_Rep1_dm6_BeS.hic_LD_DRE_Rep1_dm6_BeS_1","hic.hic_LD_DRE_Rep2_dm6_BeS.hic_LD_DRE_Rep2_dm6_BeS_1"),
	hic_LD_gypsy1    = c('hic.hic_LD_gypsy1_Rep1_dm6_BeS.hic_LD_gypsy1_Rep1_dm6_BeS_1','hic.hic_LD_gypsy1_Rep2_dm6_BeS.hic_LD_gypsy1_Rep2_dm6_BeS_1','hic.hic_LD_gypsy1_Rep3_dm6_BeS.hic_LD_gypsy1_Rep3_dm6_BeS_1'),
        hic_LD_gypsy2Enh = c('hic.hic_LD_gypsy2Enh_Rep1_dm6_BeS.hic_LD_gypsy2Enh_Rep1_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1','hic.hic_LD_gypsy2Enh_Rep3_dm6_BeS.hic_LD_gypsy2Enh_Rep3_dm6_BeS_1'),
        hic_LD_gypsy3    = c('hic.hic_LD_gypsy3_Rep1_dm6_SD.hic_LD_gypsy3_Rep1_dm6_SD_1','hic.hic_LD_gypsy3_Rep1res_dm6_SD.hic_LD_gypsy3_Rep1res_dm6_SD_1','hic.hic_LD_gypsy3_Rep2_dm6_SD.hic_LD_gypsy3_Rep2_dm6_SD_1','hic.hic_LD_gypsy3_Rep2res_dm6_SD.hic_LD_gypsy3_Rep2res_dm6_SD_1')
	)

registerDoParallel(cores=12)

sampleCoverageAll <- list()
sampleInts        <- list()
targetInts        <- list()
refTrackNames     <- list()

chr <- 'chr2L'

for(refTrack in names(refData))
{
    for(currentTrack in refData[[refTrack]])
    {
	sampleCoverage <- 0
	print(refTrack)

	refTrackName  <- gsub('_merge_mm10_IJ','',gsub('chic_','',refTrack))
	refTrackNames[[refTrack]] <- refTrackName

        splitName <- unlist(strsplit(currentTrack,'\\.'))
	setName <- splitName[length(splitName)]
	setName <- gsub('_merge_dm6_BS','',gsub('chic_','',setName))

	### Whole chromosomes ###
	for (chrom1 in gintervals.all()$chrom)
	{
            for (chrom2 in gintervals.all()$chrom)
	    {
		targetRegion2D <- gintervals.2d(chrom1,0,gintervals.all()[gintervals.all()$chrom == chrom1,'end'],chrom2,0,gintervals.all()[gintervals.all()$chrom == chrom2,'end'])

		targetData <- gextract(currentTrack,targetRegion2D,colnames=("obs"))
		if(chrom1 == chrom2)
		{
		    targetInts[[setName]] <- sum(targetData$obs)/2 # /2 because we extract both the pair i-j and the pair j-i
		} else {
		    targetInts[[setName]] <- sum(targetData$obs)/2 # /2 because we interrogate both the chrom pair chrom1-chrom2 and chrom2-chrom1
		}
		print(paste0(setName," ",chrom1," ",chrom2," ",targetInts[[setName]]))
		sampleCoverage <- sampleCoverage + targetInts[[setName]]
	    }
	}
	sampleInts[[refTrack]] <- sampleCoverage
        print(paste0(setName," ",sampleInts[[refTrack]]))

	#### Dac region ###
	targetRegion2D <- g2d(gintervals(chr,10e6,20e6))

	targetData <- gextract(currentTrack,targetRegion2D,colnames=("obs"))
	targetInts[[setName]] <- sum(targetData$obs)/2

	print(paste0(setName," dac-region ",targetInts[[setName]]))

    }
}
print(sampleInts)
