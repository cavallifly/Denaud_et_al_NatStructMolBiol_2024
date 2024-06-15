### What do you need to run the analysis and obtain the figures? (AKA Dependencies) ###
#misha
#shaman
#ggpubr
#dplyr
#doParallel
#bash: identify, convert

scriptsDir=${PWD}/scripts/

# list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed
samples="LE_WT larvae_DWT pupae_DWT LD_gypsy1 LD_gypsy2 LD_gypsy3"

# Mapping of the datasets and the mishaDB entries for the observed tracks
# have been obtained using the "shHiC2" pipeline in
# T. Nagano, Y. Lubling, C. Várnai, C. Dudley, W. Leung, Y. Baran, N. M. Cohen, S. Wingett, P. Fraser, A. Tanay,
# Cell-cycle dynamics of chromosomal organization at single-cell resolution. Nature 547, 61–67 (2017).

# Compute statistics on the mapped samples
#bash     ${scriptsDir}/01_computeValidPairsFromTracks.sh &> 01_computeValidPairsFromTracks.out
#bash     ${scriptsDir}/01_count_readPairs_from_fastq.sh  &> 01_count_readPairs_from_fastq.out

# Compute the shuffled (expected) tracks and the shaman Hi-C scores
#for sample in ${samples} ;
#do
#    Rscript ${scriptsDir}/01_generateShuffleTrack.R ${sample} &>> 01_generateShuffleTrack.out
#    Rscript ${scriptsDir}/02_computeScoreTrack.R    ${sample} &>> 02_computeScoreTrack.out
#done # Close cycle over ${sample}

# Analysis of the Hi-C maps
#Rscript ${scriptsDir}/04_computeInsulationTracks.R  &> 04_computeInsulationTracks.out
#Rscript ${scriptsDir}/05_computeDownSampledTracks.R &> 05_computeDownSampledTracks.out
#Rscript ${scriptsDir}/06_computeDiffScoreTrack_parallel.R &> 06_computeDiffScoreTrack_parallel.out
###Rscript ${scriptsDir}/06a_plotObsMaps.R &> 06a_plotObsMaps.out
###Rscript ${scriptsDir}/06b_plotICEMaps.R &> 06b_plotICEMaps.out

# Data and Figures for each figure panel
Rscript ${scriptsDir}/07_ScoreMapsPlot_Fig2a.R                 &> 07_ScoreMapsPlot_Fig2a.out                   
Rscript ${scriptsDir}/08a_get_obsContacts_PcG_domains_Fig2b.R  &>> 08a_get_obsContacts_PcG_domains_Fig2b.out   
Rscript ${scriptsDir}/08b_plot_obsContacts_PcG_domains_Fig2b.R &>> 08b_plot_obsContacts_PcG_domains_Fig2b.out  
Rscript ${scriptsDir}/09_ScoreMapsQuantification_Fig2b.R       &> 09_ScoreMapsQuantification_Fig2b.out         
Rscript ${scriptsDir}/10_ScoreMapsPlot_Fig3c.R                 &> 10_ScoreMapsPlot_Fig3c.out
Rscript ${scriptsDir}/11_ScoreMapsQuantification_Fig3d.R       &> 11_ScoreMapsQuantification_Fig3d.out
Rscript ${scriptsDir}/12_INSplot_Fig3e.R                       &> 12_INSplot_Fig3e.out
Rscript ${scriptsDir}/13_ScoreMapsPlot_Fig4a.R                 &> 13_ScoreMapsPlot_Fig4a.out
Rscript ${scriptsDir}/14_DiffScoreMapsPlot_Fig4b.R             &> 14_DiffScoreMapsPlot_Fig4b.out
Rscript ${scriptsDir}/15_INSplot_Fig4c.R                       &> 15_INSplot_Fig4c.out
Rscript ${scriptsDir}/16_INSquantification_Fig4d.R             &> 16_INSquantification_Fig4d.out
Rscript ${scriptsDir}/17_ScoreMapsQuantification_Fig4e.R       &> 17_ScoreMapsQuantification_Fig4e.out

# Obtain the triangular maps
bash ${scriptsDir}/18_generate_triangular_pngs.sh &> 18_generate_triangular_pngs.out

# Move results in folders
mkdir -p Data_for_figures
mkdir -p Figure_panels

mv *.png *.pdf Figure_panels
mv *.tab *.tsv Data_for_figures/
