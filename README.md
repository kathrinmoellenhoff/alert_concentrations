# alert_concentrations

This repository contains the code and the simulated data sets corresponding to the manuscript: "Identifying alert concentrations using a model-based bootstrap approach".

The folder "Data" contains the simulated data sets of all 18 scenarios as RData-files. They can easily be uploaded in R. 

The R script alert.R contains R code for applying the procedure as described in the manuscript (and reproducing the simulation results). The user can specify the functional form (currently available: sigEmax model and beta model).

If Linux or MacOS is used, the user can easily choose parallelization of the procedure. Therefore only the %do% command in line xx has to be changed to %dopar%. No other changes are required. But note that this parallelization does not work for Windows. 

