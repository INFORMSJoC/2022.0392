This folder contains the script files used to generate the experiment data included in Section 7.2 of the paper. Brief descriptions on the R files follow:

* MultiAssetERM_helpers.R: This file contains a collection of helper functions needed to run the simulation experiments in MultiAssetERM_main.R.
* MultiAssetERM_main.R: This script sets up the multi-asset ERM problem in Section 7.2 of the paper then runs the optimal LR nested simulation experiment design as well as SNS+, SNS, and regression-based methods. Results are saved to a .RData file.
  * This script is set up to be run (in parallel) on a linux server. See runfile.sh for a sample shell script.
  * User should set the number of macro replications (num_macro) in this script and set the number of seeds (i.e., number of parallel threads) in the shell script.
  * User should set a working directory (in line 20) where the results will be saved. Each parallel thread will produce a separate .RData file. These data files will be processed and combined by MultiAssetERM_postprocessing.R
* MultiAssetERM_postprocessing.R: This script loads and processes the simulation results produced by MultiAssetERM_main.R. It combines the different .RData files with different seeds. This script produces Table 2 and Table 3 in Section 7.2 of the paper.
