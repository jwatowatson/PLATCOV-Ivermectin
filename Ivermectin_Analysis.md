---
title: "PLATCOV Statistical Analysis: Ivermectin"
author: "James Watson"
date: "09 July, 2022"
output:
  html_document:
    toc: yes
    fig_caption: yes
    keep_md: yes
  pdf_document:
    toc: yes
  word_document:
    toc: yes
---





## Preambule

This is the analysis of the Ivermectin arm of the PLATCOV study

Data preparation is done in a different R script called data_prep.R.
This Markdown script assumes that the data are saved in a .csv file *Ivermectin_analysis.csv* in long format. 
The intention to treat population (all randomised patients included in this analysis)
The file *Ivermectin_analysis.csv* contains the patient clinical and viral load data with the following column headers:

- ID: anonymised patient id code
- Time: time from randomisation
- Trt: treatment allocation as written in CRF
- Site: site at enrolment
- Timepoint_ID: Day of study (integer 0 to 14)
- Swab_ID: RTS or TSL (right versus left tonsil)
- Plate: unique Plate ID for the PCR assay (matching with plate identifiers in interim_control_dat.csv)
- Rand_date: date of randomisation
- Any_dose: (0/1) any doses of any manufacturer received
- N_dose: integer number of doses received (any manufacturer)
- Antibody_test: 0/1 (negative/positive for SARS-CoV-2 antibody rapid test)
- Weight (kg)
- BMI: kg/weight^2
- Age: (years - has to be between 18-50)
- Sex: 0/1 (male: 1; female/other: 0)
- Symptom_onset: time since onset of symptoms (days)
- Variant: variant of concern (using standard WHO terminology for the main lineages, reference will be the predominant variant in the dataset at the start of the study)
- CT_NS: observed CT value for the N/S gene
- CT_RNaseP: observed CT value for the human RNase P gene
- Per_protocol_sample: whether at the time of sampling the patient was still in per protocol with respect to drug dosing
- log10_viral_load: log10 number of viral copies per mL (estimated from control samples using a mixed effects model)
- log10_cens_vl: censoring value



## Computational setup


```
##                _                           
## platform       x86_64-apple-darwin17.0     
## arch           x86_64                      
## os             darwin17.0                  
## system         x86_64, darwin17.0          
## status                                     
## major          4                           
## minor          0.2                         
## year           2020                        
## month          06                          
## day            22                          
## svn rev        78730                       
## language       R                           
## version.string R version 4.0.2 (2020-06-22)
## nickname       Taking Off Again
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] plotrix_3.8-2        dplyr_1.0.7          reshape2_1.4.4      
##  [4] tictoc_1.0.1         censReg_0.5-32       maxLik_1.5-2        
##  [7] miscTools_0.6-26     RColorBrewer_1.1-2   loo_2.4.1           
## [10] rstanarm_2.21.1      Rcpp_1.0.7           lme4_1.1-27.1       
## [13] Matrix_1.3-4         rstan_2.21.2         ggplot2_3.3.5       
## [16] StanHeaders_2.21.0-7
## 
## loaded via a namespace (and not attached):
##   [1] plm_2.6-0          minqa_1.2.4        colorspace_2.0-2  
##   [4] ellipsis_0.3.2     ggridges_0.5.3     rsconnect_0.8.24  
##   [7] markdown_1.1       base64enc_0.1-3    rstudioapi_0.13   
##  [10] glmmML_1.1.1       DT_0.19            fansi_0.5.0       
##  [13] codetools_0.2-18   splines_4.0.2      knitr_1.34        
##  [16] shinythemes_1.2.0  bayesplot_1.8.1    Formula_1.2-4     
##  [19] jsonlite_1.7.2     nloptr_1.2.2.2     shiny_1.6.0       
##  [22] compiler_4.0.2     backports_1.2.1    assertthat_0.2.1  
##  [25] fastmap_1.1.0      cli_3.0.1          later_1.3.0       
##  [28] htmltools_0.5.2    prettyunits_1.1.1  tools_4.0.2       
##  [31] igraph_1.2.6       gtable_0.3.0       glue_1.4.2        
##  [34] V8_3.4.2           jquerylib_0.1.4    vctrs_0.3.8       
##  [37] nlme_3.1-153       crosstalk_1.1.1    lmtest_0.9-38     
##  [40] xfun_0.26          stringr_1.4.0      rbibutils_2.2.7   
##  [43] ps_1.6.0           collapse_1.7.6     mime_0.11         
##  [46] miniUI_0.1.1.1     lifecycle_1.0.0    gtools_3.9.2      
##  [49] MASS_7.3-54        zoo_1.8-9          scales_1.1.1      
##  [52] colourpicker_1.1.0 promises_1.2.0.1   parallel_4.0.2    
##  [55] sandwich_3.0-1     inline_0.3.19      shinystan_2.5.0   
##  [58] yaml_2.2.1         curl_4.3.2         gridExtra_2.3     
##  [61] sass_0.4.0         bdsmatrix_1.3-4    stringi_1.7.4     
##  [64] dygraphs_1.1.1.6   checkmate_2.0.0    boot_1.3-28       
##  [67] pkgbuild_1.2.0     Rdpack_2.1.3       rlang_0.4.11      
##  [70] pkgconfig_2.0.3    matrixStats_0.61.0 evaluate_0.14     
##  [73] lattice_0.20-44    purrr_0.3.4        rstantools_2.1.1  
##  [76] htmlwidgets_1.5.4  processx_3.5.2     tidyselect_1.1.1  
##  [79] plyr_1.8.6         magrittr_2.0.1     R6_2.5.1          
##  [82] generics_0.1.0     DBI_1.1.1          pillar_1.6.2      
##  [85] withr_2.4.2        xts_0.12.1         survival_3.2-13   
##  [88] tibble_3.1.4       crayon_1.4.1       utf8_1.2.2        
##  [91] rmarkdown_2.11     grid_4.0.2         callr_3.7.0       
##  [94] threejs_0.3.3      digest_0.6.27      xtable_1.8-4      
##  [97] httpuv_1.6.3       RcppParallel_5.1.4 stats4_4.0.2      
## [100] munsell_0.5.0      bslib_0.3.0        shinyjs_2.0.0
```

## Load data


```
## ITT population:
```

```
## 
##   Favipiravir    Fluoxetine    Ivermectin No study drug     Regeneron 
##            45             3            46            45            40 
##    Remdesivir 
##            45
```

```
## These IDs are in the ITT database but are not in the PCR database:
##   [1] "PLT-TH1-003"  "PLT-TH1-005"  "PLT-TH1-006"  "PLT-TH1-012"  "PLT-TH1-013" 
##   [6] "PLT-TH1-016"  "PLT-TH1-019"  "PLT-TH1-022"  "PLT-TH1-025"  "PLT-TH1-027" 
##  [11] "PLT-TH1-029"  "PLT-TH1-030"  "PLT-TH1-031"  "PLT-TH1-033"  "PLT-TH1-039" 
##  [16] "PLT-TH1-042"  "PLT-TH1-044"  "PLT-TH1-046"  "PLT-TH1-048"  "PLT-TH1-052" 
##  [21] "PLT-TH1-054"  "PLT-TH1-055"  "PLT-TH1-056"  "PLT-TH1-057"  "PLT-TH1-059" 
##  [26] "PLT-TH1-060"  "PLT-TH1-061"  "PLT-TH1-062"  "PLT-TH1-063"  "PLT-TH1-064" 
##  [31] "PLT-TH1-065"  "PLT-TH1-067"  "PLT-TH1-070"  "PLT-TH1-071"  "PLT-TH1-073" 
##  [36] "PLT-TH1-074"  "PLT-TH1-076"  "PLT-TH1-078"  "PLT-TH1-079"  "PLT-TH1-082" 
##  [41] "PLT-TH1-083"  "PLT-TH1-084"  "PLT-TH1-087"  "PLT-TH1-089"  "PLT-TH1-090" 
##  [46] "PLT-TH1-091"  "PLT-TH1-092"  "PLT-TH1-095"  "PLT-TH1-096"  "PLT-TH1-097" 
##  [51] "PLT-TH1-100"  "PLT-TH1-101"  "PLT-TH1-102"  "PLT-TH1-103"  "PLT-TH1-106" 
##  [56] "PLT-TH1-108"  "PLT-TH1-109"  "PLT-TH1-110"  "PLT-TH1-112"  "PLT-TH1-113" 
##  [61] "PLT-TH1-114"  "PLT-TH1-117"  "PLT-TH1-118"  "PLT-TH1-119"  "PLT-TH1-121" 
##  [66] "PLT-TH1-122"  "PLT-TH1-125"  "PLT-TH1-126"  "PLT-TH1-127"  "PLT-TH1-130" 
##  [71] "PLT-TH1-131"  "PLT-TH1-132"  "PLT-TH1-135"  "PLT-TH1-136"  "PLT-TH1-138" 
##  [76] "PLT-TH1-140"  "PLT-TH1-143"  "PLT-TH1-145"  "PLT-TH1-146"  "PLT-TH1-147" 
##  [81] "PLT-TH1-148"  "PLT-TH1-150"  "PLT-TH1-151"  "PLT-TH1-152"  "PLT-TH1-153" 
##  [86] "PLT-TH1-154"  "PLT-TH1-155"  "PLT-TH1-163"  "PLT-TH1-164"  "PLT-TH1-166" 
##  [91] "PLT-TH1-167"  "PLT-TH1-169"  "PLT-TH1-170"  "PLT-TH1-171"  "PLT-TH1-172" 
##  [96] "PLT-TH1-174"  "PLT-TH1-177"  "PLT-TH1-178"  "PLT-TH1-181"  "PLT-TH1-182" 
## [101] "PLT-TH1-184"  "PLT-TH1-186"  "PLT-TH1-190"  "PLT-TH1-192"  "PLT-TH1-193" 
## [106] "PLT-TH1-194"  "PLT-TH1-196"  "PLT-TH1-197"  "PLT-TH1-198"  "PLT-TH1-199" 
## [111] "PLT-TH1-200"  "PLT-TH1-201"  "PLT-TH1-203"  "PLT-TH1-204"  "PLT-TH58-001"
## [116] "PLT-TH58-003" "PLT-TH58-005" "PLT-TH58-007" "PLT-TH58-009" "PLT-TH57-004"
## [121] "PLT-TH57-005" "PLT-TH57-006" "PLT-TH57-007" "PLT-TH57-008"
```

```
## [1] TRUE
```

```
## Negative time for following samples: PLT-TH57-003
## Negative time for following samples: PLT-TH57-009
## Negative time for following samples: PLT-TH57-010
## Negative time for following samples: PLT-TH58-004
```

```
## All negative samples for id: PLT-TH1-128
```


## Data summaries

Display the per protocol matrix


```
## Number of patients per arm in modified intention to treat analysis
```

```
##                          Include_mITT
##                           FALSE TRUE
##   Casirivimab/\nimdevimab     0   10
##   Ivermectin                  1   45
##   No study drug               3   41
```

```
## Number of swabs per protocol per treatment
```

```
##                          PP_swabs
##                            6  8 12 20
##   Casirivimab/\nimdevimab  0  0  0 10
##   Ivermectin               1  1  3 41
##   No study drug            2  0  0 42
```




```
## We have 1992 PCR datapoints on 99 patients from 3 sites between 2021-09-30 and 2022-04-18
```

![](Ivermectin_Analysis_files/figure-html/data_summaries-1.png)<!-- -->

```
## [1] TRUE
```

```
## The analysis dataset contains 96 patients. The geometric mean baseline (defined as samples taken within 6 hours of randomisation) viral load was 361704 copies per mL (IQR: 78075 to 2777174; range from 63 to 80411194)
```


Summary table


```
## [1] "Ivermectin"              "Casirivimab/\nimdevimab"
## [3] "No study drug"
```



Table: Summary of patient characteristics included in the current interim analysis (n= 96). Age: median (range); baseline viral load (log10 copies per mL: mean (range)); vaccinated: % with any number of doses; number of vaccine doses: median (range); antibody data are from rapid tests done at screening (+ is presence of IgM or IgG band).

|Arm                    |  n|Age          |Baseline viral load (log10) |Number of vaccine doses | Antibody+ (%)| Male (%)| th001| th057| th058|
|:----------------------|--:|:------------|:---------------------------|:-----------------------|-------------:|--------:|-----:|-----:|-----:|
|Casirivimab/
imdevimab | 10|26.5 (18-31) |5.5 (3.7-7.8)               |2 (0-3)                 |            50|       20|    10|     0|     0|
|Ivermectin             | 45|29 (19-45)   |5.7 (1.9-7.6)               |2 (0-4)                 |            78|       47|    41|     2|     2|
|No study drug          | 41|27 (20-43)   |5.5 (3-7.7)                 |2 (2-4)                 |            90|       44|    36|     3|     2|


## Summary data plot


```
## Plotting data for 96 individuals
```

![](Ivermectin_Analysis_files/figure-html/trt_data_plot-1.png)<!-- -->


## Model fitting
### Specify priors




### Prepare model

Make stan data set.

Covariates that we use in model 2:

* Vaccination (number of doses)
* Age (standardised to have mean=0 and sd=1)
* Time since symptom onset (days, between 0 and 4)
* Variant (WHO variants of concern) 
* Serology rapid test (+/-)



```
## Total number of datapoints up until day 8 is 1700
```

```
## Number of patients per arm in analysis:
```

```
## 
##           No study drug Casirivimab/\nimdevimab              Ivermectin 
##                      41                      10                      45
```

```
## There are a total of 96 patients in the database with a total of 1700 PCRs analysable
## 7.12% of samples are below LOD
## check stan data formatting:
```


### Setup model runs

We fit a set of Bayesian hierarchical models.

There are three underlying stan models
* *Linear_model_basic.stan*: vanilla student-t regression with left censoring at 0 and with individual random effects for slope and intercept;
* *Linear_model_RNaseP.stan*: Same as before but with the RNaseP measurements;
* *Nonlinear_model_RNaseP.stan*: Non-linear model (up and then down) with RNaseP adjustment.

Models 2 and 3 are combined with either informative priors or non-informative priors, and with or without full covariate adjustment (8 combinations). Model 1 is only run with informative priors and with only key covariates.


```
## [1] "Stan_models/Linear_model_basic.stan"    
## [2] "Stan_models/Linear_model_RNaseP.stan"   
## [3] "Stan_models/Nonlinear_model_RNaseP.stan"
```

```
## We are running all models with 4 chains and 5000 samples for each chain, discarding half for burn-in and thining every 10, thus giving a total of 1200 posterior samples per model.
```


Models are run on a remote cluster using the R script *run_models.R'* and the bash script *bmrc.sh*. Each model is given a seed for reproducibility.

Load model fits:



### Model fits: summaries





### Model comparisons using loo


```
## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2563.1 32.6
## p_loo       168.5  4.5
## looic      5126.2 65.2
## ------
## Monte Carlo SE of elpd_loo is 0.5.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1693  99.6%   426       
##  (0.5, 0.7]   (ok)          7   0.4%   704       
##    (0.7, 1]   (bad)         0   0.0%   <NA>      
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## 
## All Pareto k estimates are ok (k < 0.7).
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2537.0 32.2
## p_loo       169.6  4.6
## looic      5074.1 64.3
## ------
## Monte Carlo SE of elpd_loo is 0.5.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1698  99.9%   404       
##  (0.5, 0.7]   (ok)          2   0.1%   488       
##    (0.7, 1]   (bad)         0   0.0%   <NA>      
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## 
## All Pareto k estimates are ok (k < 0.7).
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2509.7 32.6
## p_loo       216.3  6.1
## looic      5019.4 65.1
## ------
## Monte Carlo SE of elpd_loo is NA.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1676  98.6%   219       
##  (0.5, 0.7]   (ok)         23   1.4%   215       
##    (0.7, 1]   (bad)         1   0.1%   319       
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2537.5 32.2
## p_loo       169.5  4.5
## looic      5074.9 64.4
## ------
## Monte Carlo SE of elpd_loo is 0.5.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1697  99.8%   445       
##  (0.5, 0.7]   (ok)          3   0.2%   431       
##    (0.7, 1]   (bad)         0   0.0%   <NA>      
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## 
## All Pareto k estimates are ok (k < 0.7).
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2510.2 32.7
## p_loo       217.5  6.1
## looic      5020.4 65.3
## ------
## Monte Carlo SE of elpd_loo is NA.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1684  99.1%   261       
##  (0.5, 0.7]   (ok)         11   0.6%   130       
##    (0.7, 1]   (bad)         4   0.2%   141       
##    (1, Inf)   (very bad)    1   0.1%   288       
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2537.3 32.3
## p_loo       170.6  4.6
## looic      5074.7 64.5
## ------
## Monte Carlo SE of elpd_loo is 0.5.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1698  99.9%   374       
##  (0.5, 0.7]   (ok)          2   0.1%   1156      
##    (0.7, 1]   (bad)         0   0.0%   <NA>      
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## 
## All Pareto k estimates are ok (k < 0.7).
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2510.7 32.7
## p_loo       218.3  6.1
## looic      5021.3 65.4
## ------
## Monte Carlo SE of elpd_loo is NA.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1675  98.5%   352       
##  (0.5, 0.7]   (ok)         21   1.2%   132       
##    (0.7, 1]   (bad)         4   0.2%   338       
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2538.9 32.3
## p_loo       172.1  4.6
## looic      5077.8 64.7
## ------
## Monte Carlo SE of elpd_loo is 0.5.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1696  99.8%   332       
##  (0.5, 0.7]   (ok)          4   0.2%   769       
##    (0.7, 1]   (bad)         0   0.0%   <NA>      
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## 
## All Pareto k estimates are ok (k < 0.7).
## See help('pareto-k-diagnostic') for details.
```

```
## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```
## 
## Computed from 1200 by 1700 log-likelihood matrix
## 
##          Estimate   SE
## elpd_loo  -2512.8 32.8
## p_loo       220.7  6.2
## looic      5025.7 65.6
## ------
## Monte Carlo SE of elpd_loo is NA.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     1675  98.5%   204       
##  (0.5, 0.7]   (ok)         24   1.4%   134       
##    (0.7, 1]   (bad)         1   0.1%   936       
##    (1, Inf)   (very bad)    0   0.0%   <NA>      
## See help('pareto-k-diagnostic') for details.
```

```
##        elpd_diff se_diff
## model2   0.0       0.0  
## model1 -27.4       9.0
```

```
##        elpd_diff se_diff
## model2   0.0       0.0  
## model1 -27.2       9.0
```

```
##        elpd_diff se_diff
## model2   0.0       0.0  
## model1 -26.7       9.1
```

```
##        elpd_diff se_diff
## model2   0.0       0.0  
## model1 -26.1       9.0
```


## Results

### Estimated treatment effects under the 4 models

Posterior distributions over the treatment effects for the interventions. Red: no effect; blue: median inferred effect.


```
## 
## *******************
## Mean estimated treatment effects (multiplicative):
```

```
##       Casirivimab/\nimdevimab Ivermectin
##  [1,]                1.548114  0.9034610
##  [2,]                1.523211  0.9085989
##  [3,]                1.456219  0.9079275
##  [4,]                1.655142  0.9216786
##  [5,]                1.529204  0.9089353
##  [6,]                1.511075  0.9198498
##  [7,]                1.442625  0.9107858
##  [8,]                1.649500  0.9306607
##  [9,]                1.494531  0.9115989
```

```
## 
## *******************
## Probability of super-superiority:
```

```
##       Casirivimab/\nimdevimab Ivermectin
##  [1,]                    94.8        2.7
##  [2,]                    95.5        2.1
##  [3,]                    93.2        1.5
##  [4,]                    97.5        4.2
##  [5,]                    95.1        1.9
##  [6,]                    92.9        2.7
##  [7,]                    92.2        1.6
##  [8,]                    97.2        4.6
##  [9,]                    92.7        1.9
```


Overall effects


```
##                           mean         sd        2.5%       97.5%     n_eff
## trt_effect[1]       0.42082037 0.18000138  0.06723114  0.76607278 1215.0412
## trt_effect[2]      -0.09585150 0.11045454 -0.31760287  0.11183423  714.9721
## alpha_0             5.70475575 0.16384894  5.37658822  6.00774713 1072.5793
## beta_0             -0.38285654 0.04776558 -0.48664422 -0.30260132 1060.9486
## sigma_logvl         0.92111378 0.02860333  0.86297974  0.97579958 1203.0559
## sigmasq_u[1]        0.94172296 0.08192802  0.78784124  1.11467205 1252.0824
## sigmasq_u[2]        0.49375166 0.05594806  0.39692484  0.61615426 1244.5064
## t_dof               6.36524522 1.03705028  4.70147254  8.68303784 1218.9419
## gamma_rnasep        0.21260111 0.02898213  0.15517090  0.26782564 1264.5385
## slope_coefs[1]      0.04646932 0.15304557 -0.24656406  0.35617377 1259.6781
## slope_coefs[2]      0.03347912 0.13870962 -0.24348413  0.28651548 1040.3141
## slope_coefs[3]      0.34013473 0.21968844 -0.09551169  0.78132210 1266.0676
## slope_coefs[4]      0.06459665 0.25177572 -0.44811509  0.54125861 1227.7465
## intercept_coefs[1] -0.44623992 0.24176494 -0.92683279  0.03669462 1163.8724
## intercept_coefs[2]  0.55373702 0.21690775  0.13535158  0.97477534 1083.4873
## intercept_coefs[3] -0.28655870 0.33925957 -0.97229145  0.41115809 1315.7318
## intercept_coefs[4] -0.34183028 0.37350739 -1.03969406  0.38853144 1237.8437
##                         Rhat
## trt_effect[1]      0.9987007
## trt_effect[2]      1.0056810
## alpha_0            1.0022903
## beta_0             1.0013694
## sigma_logvl        1.0019461
## sigmasq_u[1]       0.9992692
## sigmasq_u[2]       0.9990045
## t_dof              1.0003107
## gamma_rnasep       0.9979734
## slope_coefs[1]     1.0014472
## slope_coefs[2]     1.0033287
## slope_coefs[3]     1.0002600
## slope_coefs[4]     1.0007402
## intercept_coefs[1] 1.0005099
## intercept_coefs[2] 1.0004774
## intercept_coefs[3] 0.9989030
## intercept_coefs[4] 1.0005374
```

```
## Under model 2 the change in rate of clearance for Ivermectin compared to no study drug is -9.14% (95%CI: -27.21% to 11.83%)
```

```
## Under model 2 the change in rate of clearance for Regeneron compared to no study drug is 52.32% (95%CI: 6.95% to 115.13%)
```

```
## The main 3 models with weakly informative priors:
```

![](Ivermectin_Analysis_files/figure-html/treatment_effect_plot-1.png)<!-- -->


Covariate effects on the intercept (baseline viral load) and slope (viral clearance):


```
## The following `from` values were not present in `x`: Age_scaled, Antibody_test, Symptom_onset, N_dose
## The following `from` values were not present in `x`: Age_scaled, Antibody_test, Symptom_onset, N_dose
```

![](Ivermectin_Analysis_files/figure-html/cov_effects-1.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/cov_effects-2.png)<!-- -->

### Slopes over time

Plot the absolute slope estimate for each individual over time


```
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
```

```
## Slopes plot for model setting 2
```

```
##                                    mod prior cov_matrices Niter Nwarmup Nthin
## 2 Stan_models/Linear_model_RNaseP.stan     1            1  5000    2000    10
##   Nchain
## 2      4
```

```
## In the no study drug arm the mean clearance half life was 20.1 (range 7.1 to 41.4)
```

```
## In the Ivermectin arm the mean clearance half life was 21.5 (range 8.9 to 51.9)
```

```
## In the Regeneron arm the mean clearance half life was 12.7 (range 6.2 to 19.9)
```

```
## The model estimated population mean clearance half-life is 19.2 (95% CI 14.8-23.9)
```

![](Ivermectin_Analysis_files/figure-html/slopes_over_time-1.png)<!-- -->


changes in half life


```
## In the Ivermectin arm the mean change in half life is 1.93 (95% CI -2.13 to 6.57)
```

```
## In the Regeneron arm the mean change in half life is -6.53 (95% CI -12.01 to -1.06)
```


Plot the individual slope estimates by group


```
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
```

![](Ivermectin_Analysis_files/figure-html/slopes_by_group-1.png)<!-- -->


Some exploratory covariate analyses


```
## 
## Call:
## lm(formula = t_12 ~ trt, data = trt_summary_dat)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -13.071  -5.065  -1.746   4.714  30.331 
## 
## Coefficients:
##                            Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  20.138      1.289  15.621   <2e-16 ***
## trtCasirivimab/\nimdevimab   -7.412      2.911  -2.546   0.0125 *  
## trtIvermectin                 1.399      1.782   0.785   0.4344    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.255 on 93 degrees of freedom
## Multiple R-squared:  0.09123,	Adjusted R-squared:  0.07169 
## F-statistic: 4.668 on 2 and 93 DF,  p-value: 0.0117
```

![](Ivermectin_Analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```
## 
## Call:
## lm(formula = rate_mean ~ AUC_ivm, data = trt_summary_dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.61810 -0.06836  0.02267  0.12154  0.26675 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -4.042e-01  2.550e-02 -15.848   <2e-16 ***
## AUC_ivm     -1.910e-07  3.298e-06  -0.058    0.954    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1694 on 84 degrees of freedom
##   (10 observations deleted due to missingness)
## Multiple R-squared:  3.991e-05,	Adjusted R-squared:  -0.01186 
## F-statistic: 0.003353 on 1 and 84 DF,  p-value: 0.954
```

```
## 
## Call:
## lm(formula = rate_mean ~ Cmax_ivm, data = trt_summary_dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.61856 -0.06838  0.02325  0.12166  0.26689 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -4.037e-01  2.551e-02 -15.824   <2e-16 ***
## Cmax_ivm    -9.336e-06  1.117e-04  -0.084    0.934    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1694 on 84 degrees of freedom
##   (10 observations deleted due to missingness)
## Multiple R-squared:  8.314e-05,	Adjusted R-squared:  -0.01182 
## F-statistic: 0.006984 on 1 and 84 DF,  p-value: 0.9336
```


Illustrative PK plot

![](Ivermectin_Analysis_files/figure-html/pk_predicted_profiles-1.png)<!-- -->


### Individual plots

Individual plots colored by model


```
## The following `from` values were not present in `x`: PLT-TH1-018, PLT-TH1-020, PLT-TH1-105
```

![](Ivermectin_Analysis_files/figure-html/individ_fits-1.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-2.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-3.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-4.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-5.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-6.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/individ_fits-7.png)<!-- -->


## Sensitivity analysis

### Treatment effects

![](Ivermectin_Analysis_files/figure-html/treatment_effect_sensitivity-1.png)<!-- -->


### Left vs right tonsil

![](Ivermectin_Analysis_files/figure-html/left_versus_right-1.png)<!-- -->![](Ivermectin_Analysis_files/figure-html/left_versus_right-2.png)<!-- -->


