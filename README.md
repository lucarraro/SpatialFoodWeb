# SpatialFoodWeb

R code supporting "Coupled biological and hydrological processes shape spatial food-web structures in riverine metacommunities", by Hsi-Cheng Ho, Florian Altermatt, Luca Carraro.  

## Content

Files are described following the logical order under which they should be run.

- `create_FW_OCN.R`: creates 100 meta-food webs (saved in `utilities/100FW_nSp100.mat`) and the OCN (saved in `utilities/OCN.rda` and `utilities/OCN.mat`).
- `RUN_SFW.m`: executes the SFW model for default parameter values and performs the sensitivity analysis. Results are saved in `utilities/results_for_R.mat`.
- `RUN_UWB.m`: executes the UWB model. Results are saved in `utilities/resultsUWB_for_R.mat`.
- `ANALYZE_DATA.R`: analyzes output from SFW and UWB models, runs RMW and RND models, and produces all manuscript figures (except Figs. 1a, 1b, 6).
- `create_Fig6.mat`: produces Fig. 6 of the manuscript.
- `results`: folder storing output from SFW model. Each file is the simulation output for a given meta-food-web realization and parameter set.
- `results_UWB`: folder storing output from UWB model. Each file is the simulation output for a given meta-food-web realization.
- `utilities`: folder storing intermediate data files produced by the scripts above, as well as functions called by the above scripts.
