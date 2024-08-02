This repository contains the scripts to reproduce the figures on the CnR analysis in the manuscript Denaud et al 2024

# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

Ensure that you have a working version of R. The scripts in this repository have been tested for version 4.0.5 (2021-03-31).
Install the following R packages: misha, shaman, ggpubr, dplyr, doParallel.
Install the following bash commands: convert (part of imagemagick) and identify.
For example, you can use the commands:
```
conda install -y bioconda::r-misha
conda install -y bioconda::r-shaman
conda install -y conda-forge::r-ggpubr
conda install -y conda-forge::r-dplyr
conda install -y conda-forge::r-doparallel
conda install -y conda-forge::imagemagick
conda install -y conda-forge::identify
```

# Input data #
Next, you should download the misha tracks from GEO for the observed counts in the directory HiC_analysis/mishaDB/trackdb/dm6/tracks/hic/ of this repository.
- GSM7888066	Larv WT HiC Repli1
- GSM7888067	Larv WT HiC Repli2
- GSM7888068	Larv WT HiC Repli3
- GSM7888069	pupe WT Repl1
- GSM7888070	pupe WT Repl2
- GSM7888071	pupe WT Repl3
- GSM7888072	Larv RE deletion HiC Repli1
- GSM7888073	Larv RE deletion HiC Repli2
- GSM7888074	Larv gypsy1 insertion HiC Repli1
- GSM7888075	Larv gypsy1 insertion HiC Repli2
- GSM7888076	Larv gypsy1 insertion HiC Repli3
- GSM7888077	Larv gypsy2 insertion HiC Repli1
- GSM7888078	Larv gypsy2 insertion HiC Repli2
- GSM7888079	Larv gypsy2 insertion HiC Repli3
- GSM7888080	Larv gypsy3 insertion HiC Repli1
- GSM7888081	Larv gypsy3 insertion HiC Repli2
- GSM7888082	Larv gypsy3 insertion HiC Repli1res
- GSM7888083	Larv gypsy3 insertion HiC Repli2res

These tracks have been obtained using the "scHiC2" [pipeline](https://github.com/tanaylab/schic2).

# Data analysis #
Now, you are ready to run the scripts in HiC_analysis/scripts. To do so, you can access the directory HiC_analysis of this repository using
```
cd HiC_analysis
```
and run one after the other the following commands:
```
scriptsDir=./scripts/

# Compute the shuffled (expected) tracks and the shaman Hi-C scores
Rscript ${scriptsDir}/01_generateShuffleTrack.R &> 02_generateShuffleTrack.out
Rscript ${scriptsDir}/02_computeScoreTrack.R    &> 03_computeScoreTrack.out

# Analysis of the Hi-C maps
Rscript ${scriptsDir}/04_computeInsulationTracks.R        &> 04_computeInsulationTracks.out
Rscript ${scriptsDir}/05_computeDownSampledTracks.R       &> 05_computeDownSampledTracks.out
Rscript ${scriptsDir}/06_computeDiffScoreTrack_parallel.R &> 06_computeDiffScoreTrack_parallel.out

# Generate data and Figures for each figure panel
Rscript ${scriptsDir}/07_ScoreMapsPlot_Fig2a.R                 &> 07_ScoreMapsPlot_Fig2a.out                   
Rscript ${scriptsDir}/08a_get_obsContacts_PcG_domains_Fig2b.R  &> 08a_get_obsContacts_PcG_domains_Fig2b.out   
Rscript ${scriptsDir}/08b_plot_obsContacts_PcG_domains_Fig2b.R &> 08b_plot_obsContacts_PcG_domains_Fig2b.out  
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

# Move results in the data and figures folders
mkdir -p Data_for_figures
mkdir -p Figure_panels

mv *.png *.pdf Figure_panels
mv *.tab *.tsv Data_for_figures/
```

Once the scripts are finished, you will obtain the panels and the data points used to obtain them for all the Figures in Denaud et al NatStructMolBiol 2024.
To obtain the final version of the figures the panels have been assembled using the PowerPoint program.

## Contributions ##
The code in this repository has been developed at the [Cavalli Lab](https://www.igh.cnrs.fr/en/research/departments/genome-dynamics/chromatin-and-cell-biology) with the contributions of Giorgio L. Papadopoulos, Marco Di Stefano, and Gonzalo Sabaris.

