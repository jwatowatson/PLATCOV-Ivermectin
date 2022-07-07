# script for running all the stan models with all settings on the BMRC cluster

args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file


## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)

load('Rout/model_run_setup.RData')

Max_job = nrow(model_settings)
if(i > Max_job) stop('no model setting corresponding to job ID')

writeLines('Doing the following job:')
print(model_settings[i, ])

options(mc.cores = model_settings$Nchain[i])
stopifnot(model_settings$Nchain[i]>getDoParWorkers()) # check worker number assigned

mod = stan_model(file = as.character(model_settings$mod[i])) # compile stan model

analysis_data_stan = stan_inputs$analysis_data_stan
analysis_data_stan$trt_mat = stan_inputs$Trt_matrix
analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)

# 0 means no covariates at all
if(model_settings$cov_matrices[i]==0){
  analysis_data_stan$K_cov_intercept=0
  analysis_data_stan$x_intercept = array(dim = c(analysis_data_stan$Ntot,0))
  analysis_data_stan$K_cov_slope=0
  analysis_data_stan$x_slope = array(dim = c(analysis_data_stan$Ntot,0))
} else {
  analysis_data_stan$x_intercept = stan_inputs$cov_matrices$X_int[[model_settings$cov_matrices[i]]]
  analysis_data_stan$K_cov_intercept=ncol(analysis_data_stan$x_intercept)
  analysis_data_stan$x_slope = stan_inputs$cov_matrices$X_slope[[model_settings$cov_matrices[i]]]
  analysis_data_stan$K_cov_slope=ncol(analysis_data_stan$x_slope)
}

# sample posterior
out = sampling(mod, 
               data=c(analysis_data_stan,
                      all_priors[[model_settings$prior[i]]]),
               iter=model_settings$Niter[i],
               chain=model_settings$Nchain[i],
               thin=model_settings$Nthin[i],
               warmup=model_settings$Nwarmup[i],
               seed=i, # for reproducibility
               save_warmup = FALSE,
               pars=c('L_Omega','theta_rand_id'), # we don't save these as it takes up vast memory!
               include=FALSE)


save(out, file = paste0('Rout/model_fits_',i,'.RData'))# save output

writeLines('Finished job')

