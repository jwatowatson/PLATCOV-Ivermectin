# PLATCOV-Ivermectin

Statistical analysis of the Ivermectin arm in the PLATCOV trial.
The PLATCOV trial is registered at clinicaltrials.gov number [NCT05041907](https://clinicaltrials.gov/ct2/show/NCT05041907)

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



## Overview

This github repo provides the data and code for the statistical analysis of the Ivermectin arm in the PLATCOV trial published in [Elife](https://elifesciences.org/articles/83201). The primary analysis of the trial consists of fitting linear and non-linear hierarchical Bayesian regression models to the serial viral load measurements over time. The regression models use left-censoring for viral loads below the lower limit of detection (i.e. a Ct value of 40 or above). The viral load is expressed as the log base 10 number of copies per mL. The statistical analysis plan used at the time of analysis is given in the file *PLATCOV_SAP_v2.1_13052022.pdf*.

All models are fit to data using Monte Carlo approximation of the posterior distributions with the open access software *stan* (interface to R with *rstan*). The folder *Stan_models* contains the stan code for the three models used:

* Linear_model_basic.stan: the most basic model (no adjustment for human RNaseP)
* Linear_model_RNaseP.stan: additional adjustment for RNaseP
* Nonlinear_model_RNaseP.stan: non-linear model allowing increases in viral load at the start

The data are given in *Ivermectin_analysis.csv*.

The RMarkdown file *Ivermectin_Analysis.Rmd* does the following:

* Loads the data and outputs summary statistics and summary plots
* Sets up the parameters for all the model runs (9 in total)
* Loads the model fits and plots them

The models are run using the R script *run_models.R*. I did this on a server as it takes a while. The data dictionary is in the main RMarkdown file.

## Software needed

The R packages needed are:

* *rstan* (interfaces with *stan*)
* *loo* (for model comparison)
* *RColorBrewer* 

Any questions or comments or if any bugs spotted drop me a message at jwatowatson at gmail dot com


