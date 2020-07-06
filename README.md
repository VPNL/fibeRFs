### Differential spatial processing in ventral and lateral face-selective regions is scaffolded by structural connections
**Dawn Finzi, Jesse Gomez, Marisa Nordt, Alex A. Rezai, Sonia Poltoratski, Kalanit Grill-Spector**

This repository contains code to analyze the data, compute statistics, and make the individual figure elements. 

#### Overall organization 
Code relating to figures, statistical tests, and data analyses is provided in subdirectories `code/`, `R_code/`, and `results/`. The `code/` subdirectories contains the matlab code used to analyze the preprocessed fMRI and dMRI data and generate .mat files included in `results/`. Code for generating figures 2A (`code/study1/pRFs/s1_plotCentersBilat.m`), 4A (`code/study1/pRFs/plotAverageCoverage.m`) and 4B (`code/study1/pRFs/s1_covByRadialDist_contraOnly.m`) are also located here. Code to generate all other figures, and code for all the statistical analyses reported in the paper, are included in `R_code/`, labeled by the corresponding figure number. `results/` includes the processed data used for figure generation and statistics.
Expected outputs are .mat, .fig, .png and .pdf files. 

#### Dependencies :package:
The code in `R_code/` should be runnable on your local machine assuming the following dependencies:
```
R (>= 3.3.0)
R.matlab
ggthemes
tidyverse
plyr
lmerTest
lsmeans
dplyr
PairedData
```

The matlab analysis code in `code/` was run using Matlab2015a, and relies heavily on [mrVista](http://github.com/vistalab) and the [toonotopy processing code](https://github.com/VPNL/toonotopy). It has been tested on Mac and Linux OS only and may require changing paths to reflect your local environment. The repository should be fast to download and all scripts should run in less than a minute, unless otherwise specified within the script. Other than the specified figure generation code, the matlab code unfortunately can not be run on your local machine as it depends on absolute paths to niftis and FreeSurfer surfaces. 

Please see the data availability statement in Finzi et. al., 2020 for more and contact me with any questions at: <dfinzi@stanford.edu>
