# JohnsongrassBranchingProcess

This repository contains all the figures, the data from which plots were generated and the data generating code for the manuscript *Assessing the probability and time of weed control failure due to herbicide resistance evolution*.

## Figures and Data

For each figure in the paper there is a separate folder named **Figure_** containing up to three subfolders with all data and code necessary to recreate them:

- **Data:** data files in txt format 
- **Simulations:** Matlab functions and scripts used to generate the data
- **Analysis:** Mathematica notebook performing the analysis and generating the plots

The additional folder **StandingVariation** contains Matlab functions and scripts calculating the population composition in mutation-selection balance used as standing genetic variation for all figures except for Figure 5. 

The simulation was written in `MATLAB_R2023a` and the analysis and figure generation was implemented in `Mathematica 13.1`.
