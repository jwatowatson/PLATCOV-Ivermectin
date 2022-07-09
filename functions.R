# list of functions for plotting data

plot_effect_estimates = function(effect_ests, #list of stan outputs
                                 plot_models, # indices of models to plot in list
                                 my_pch=1,
                                 mod_cols = NULL,
                                 study_threshold){
  
  if (length(my_pch)==1) my_pch = (1:length(plot_models))+15
  if (length(my_pch)!=length(plot_models)) stop('length of my_pch needs to be the same as the number of input models')
  if (is.null(mod_cols)){
    mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
    if(length(plot_models)>9) writeLines('too many models - supply user defined colors')
  }
  if(length(mod_cols)!=length(plot_models)) stop('number of colors needs to be equal to number of models to be plotted')
  
  
  K_treatments = nrow(effect_ests[[plot_models[1]]])
  
  xlims = (exp(range(sapply(effect_ests[plot_models], rbind)))-1)*100
  x_points = pretty((xlims),6)
  plot(NA, NA, xlim = range(x_points),
       ylim = c(0.75,K_treatments+.25),
       panel.first=grid(), ylab='', yaxt='n', type='n',
       xlab = 'Change in rate of clearance (%)',
       xaxt = 'n')
  axis(1, at = x_points)
  axis(2, at = 1:K_treatments,
       labels = rownames(effect_ests[[plot_models[1]]]),
       tick = F, cex.lab=1.5, cex.axis=1.5)
  index_p = rev(seq(-.2,.2, length.out = length(plot_models)))
  abline(v=0,lwd=2)
  polygon(c(-1000, 100*(study_threshold-1), 100*(study_threshold-1), -1000),
          c(-100, -100, 100, 100), border = NA,
          col = adjustcolor('grey',.4))
  for(i in 1:length(plot_models)){
    points((exp(effect_ests[[plot_models[i]]][,'mean'])-1)*100,
           1:K_treatments+index_p[i],pch=my_pch[i],
           col=mod_cols[i],cex=1.5)
    for(j in 1:K_treatments){
      lines((exp(c(effect_ests[[plot_models[i]]][j,'2.5%'],
                   effect_ests[[plot_models[i]]][j,'97.5%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1)
      lines((exp(c(effect_ests[[plot_models[i]]][j,'10%'],
                   effect_ests[[plot_models[i]]][j,'90%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3)
    }
  }
}



plot_baseline_data = function(input_data){
  
  baseline_ind = input_data$Timepoint_ID==0
  bvl = aggregate(log10_viral_load ~ ID, input_data[baseline_ind, ], median)
  
  hist(bvl$log10_viral_load,
       breaks = seq(1,8.5,by=.5),
       xlab='Baseline viral load (RNA copies per mL)',
       ylab ='Number of patients',xlim=c(1,8.5),
       main='', xaxt ='n')
  axis(1, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  grid(); 
  hist(bvl$log10_viral_load,breaks = seq(1,8.5,by=.5),add=T)
}

plot_serial_data = function(xx,Trt, trt_cols){
  
  xx$Trt = Trt
  xx$Trt_number = as.numeric(as.factor(xx$Trt))
  
  PCR_dat = aggregate(log10_viral_load ~ ID + Timepoint_ID +
                        Trt_number + Trt, 
                      data = xx, mean)
  trt_smmry = aggregate(formula = log10_viral_load ~ Timepoint_ID+Trt_number, 
                        data = PCR_dat, FUN = median)
  PCR_dat$Timepoint_ID=jitter(PCR_dat$Timepoint_ID)
  
  gap.plot(PCR_dat$Timepoint_ID, PCR_dat$log10_viral_load,
           ylab = 'RNA copies per mL', panel.first=grid(),
           xlab = 'Time since randomization (days)',
           gap = c(7.5,13.5), gap.axis = 'x',
           yticlab = '',ytics = 2, xtics = c(0,3,6,14),
           xlim = c(0,14), type='n', yaxt='n',
           col = trt_cols[PCR_dat$Trt_number],
           ylim = c(1, max(PCR_dat$log10_viral_load)))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  IDs = unique(xx$ID)
  writeLines(sprintf('Plotting data for %s individuals', length(IDs)))
  for(id in IDs){
    ind = PCR_dat$ID==id
    gap.plot(PCR_dat$Timepoint_ID[ind],
             PCR_dat$log10_viral_load[ind],
             gap = c(7.5,13.5), gap.axis = 'x',add = T,
             col=adjustcolor(trt_cols[PCR_dat$Trt_number[ind]],.3))
  }
  gap.plot(trt_smmry$Timepoint_ID, trt_smmry$log10_viral_load,
           col= trt_cols[trt_smmry$Trt_number],
           gap = c(7.5,13.5), gap.axis = 'x',add = T,
           pch = 14+trt_smmry$Trt_number, cex=1.5)
  for(tt in unique(trt_smmry$Trt_number)){
    ind = trt_smmry$Trt_number==tt
    gap.plot(trt_smmry$Timepoint_ID[ind],
             trt_smmry$log10_viral_load[ind], 
             gap = c(7.5,13.5), gap.axis = 'x',add = T,
             type='l',col = trt_cols[tt],lwd=3)
  }
  
  trt_nums = sort(unique(PCR_dat$Trt_number))
  trt_cols = trt_cols[trt_nums]
  trts = plyr::mapvalues(trt_nums,
                         from=PCR_dat$Trt_number[!duplicated(PCR_dat$Trt_number)],
                         to=PCR_dat$Trt[!duplicated(PCR_dat$Trt_number)])
  for(i in 1:length(trts)){
    trts[i] = paste0(trts[i],' (n=',sum(PCR_dat$Trt[!duplicated(PCR_dat$ID)]==trts[i]),')')
  }
  legend('topright', col=trt_cols, 
         legend = trts, title = 'Median',
         lwd=2,pch=14+trt_nums,cex=1, inset = 0.03)
}

bayes_R2 = function(mod_preds, mod_residuals) {
  var_pred = apply(mod_preds, 1, var)
  var_res = apply(mod_residuals, 1, var)
  var_pred / (var_pred + var_res)
}


make_stan_inputs = function(input_data_fit, 
                            int_covs_base,
                            int_covs_full,
                            slope_covs_base,
                            slope_covs_full,
                            trt_frmla,
                            Dmax
){
  
  ## check censored values come last
  if(!all(diff(input_data_fit$log10_viral_load == input_data_fit$log10_cens_vl)>=0)) stop()
  ind_dup = !duplicated(input_data_fit$ID)
  
  input_data_fit$RnaseP_scaled = t(scale(40 - input_data_fit$CT_RNaseP, 
                                         scale = F))[1,] 
  input_data_fit$Age_scaled = (input_data_fit$Age-mean(input_data_fit$Age[ind_dup]))/sd(input_data_fit$Age[ind_dup])
  
  # make the covariate matrix
  # check no missing data
  if(!all(!apply(input_data_fit[, union(int_covs_full,slope_covs_full)], 2, function(x) any(is.na(x))))){
    stop('Missing data in covariate matrix!')
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_base,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_1 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_2 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_base,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_1 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_2 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  if(!nrow(X_intcpt_1) == nrow(input_data_fit)) stop()
  if(!nrow(X_slope_1) == nrow(input_data_fit)) stop()
  
  cov_matrices = list(X_int=list(X_intcpt_1, X_intcpt_2),
                      X_slope=list(X_slope_1, X_slope_2))
  
  ID_map = data.frame(ID_key = input_data_fit$ID,
                      ID_stan = as.numeric(as.factor(input_data_fit$ID)))
  writeLines(sprintf('There are a total of %s patients in the database with a total of %s PCRs analysable',
                     max(ID_map$ID_stan),
                     nrow(input_data_fit)))
  
  ind_cens = !input_data_fit$log10_viral_load>
    input_data_fit$log10_cens_vl
  
  writeLines(sprintf('%s%% of samples are below LOD',
                     round(100*mean(ind_cens),digits = 2)))
  
  analysis_data_stan = list(Ntot = nrow(input_data_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = input_data_fit$Time,
                            log_10_vl = input_data_fit$log10_viral_load,
                            log10_cens_vl = input_data_fit$log10_cens_vl,
                            RNaseP = input_data_fit$RnaseP_scaled,
                            Time_max = Dmax)
  ID_map = ID_map[!duplicated(ID_map$ID_key), ]
  
  writeLines('check stan data formatting:')
  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])
  
  Trt_matrix = model.matrix(trt_frmla, data = input_data_fit)
  Trt_matrix[,1]=0 # first column is dummy
  
  analysis_inputs = list(cov_matrices=cov_matrices,
                         analysis_data_stan=analysis_data_stan,
                         Trt_matrix=Trt_matrix,
                         ID_map=ID_map)
  return(analysis_inputs)
}


make_baseline_table = function(input_data){
  
  ind_dup = !duplicated(input_data$ID)
  platcov_sm_dat = input_data[ind_dup, ]
  
  xx_vl = aggregate(log10_viral_load ~ ID + Timepoint_ID + Trt_code, input_data, mean)
  xx_vl = aggregate(log10_viral_load ~ Trt_code, xx_vl[xx_vl$Timepoint_ID==0, ],
                    function(x) {
                      paste(round(mean(x),1),' (',
                            round(min(x),1),'-',round(max(x),1),
                            ')',sep='')})
  xx_age = aggregate(Age ~ Trt_code, platcov_sm_dat,function(x) {
    paste(round(median(x),1),' (',
          round(min(x),1),'-',round(max(x),1),
          ')',sep='')})
  xx_n = aggregate(Age ~ Trt_code, platcov_sm_dat, length)
  xx_Nvac = aggregate(N_dose ~ Trt_code, platcov_sm_dat,function(x) {
    paste(round(median(x),1),' (',
          round(min(x),1),'-',round(max(x),1),
          ')',sep='')})
  platcov_sm_dat$Vaccinated = as.numeric(platcov_sm_dat$Any_dose=='Yes')
  xx_Any_vac = aggregate(Vaccinated ~ Trt_code, platcov_sm_dat,function(x) round(100*mean(x)))
  xx_Antibody = aggregate(Antibody_test ~ Trt_code, platcov_sm_dat,function(x) round(100*mean(x)))
  xx_sex = aggregate(Sex ~ Trt_code, platcov_sm_dat,function(x) round(100*mean(x)))
  
  # Drugs allocated by Site
  xx_site = merge(data.frame(Trt_code=levels(platcov_sm_dat$Trt_code)),
                  aggregate(Site ~ Trt_code, FUN = length,data = platcov_sm_dat,
                            subset = platcov_sm_dat$Site==unique(platcov_sm_dat$Site)[1]),
                  all.x=TRUE)
  if(length(unique(platcov_sm_dat$Site)>1)){
    for(ss in unique(platcov_sm_dat$Site)[-1]){
      xx_site = cbind(xx_site,
                      merge(data.frame(Trt_code=levels(platcov_sm_dat$Trt_code)),
                            aggregate(Site ~ Trt_code, FUN = length,data = platcov_sm_dat,
                                      subset = platcov_sm_dat$Site==ss),
                            all.x=TRUE)[,-1])
    }
  }
  xx_site[is.na(xx_site)]=0
  
  xx = merge(merge(merge(merge(merge(merge(
    xx_n, xx_age,by = 'Trt_code'),
    xx_vl,by = 'Trt_code'),
    xx_Nvac,by = 'Trt_code'),
    xx_Antibody,by = 'Trt_code'),
    xx_sex, by = 'Trt_code'),
    xx_site, by = 'Trt_code')
  
  colnames(xx)=c('Arm','n','Age','Baseline viral load (log10)',
                 'Number of vaccine doses',
                 'Antibody+ (%)','Male (%)',
                 unique(platcov_sm_dat$Site))
  
  
  return(xx)
}



plot_individ_data = function(mod_out, # model fits
                             models_plot, # which models to plot
                             K_plots,
                             mod_cols,
                             ID_map,
                             analysis_data_stan
){
  
  # extract posterior parameters and outputs
  thetas = list()
  for(mm in 1:length(mod_out)){
    thetas[[mm]] = extract(mod_out[[mm]])
  }
  id = counter = 1
  
  ID_map$Trt = gsub(pattern = '\n',
                    replacement = '',
                    x = ID_map$Trt,fixed = T)
  while(id <= max(ID_map$ID_stan)){
    
    # every K_plots put a legend in bottom right panel
    # if(counter %% K_plots == 0){
    #   plot(NA,NA,xlab='',ylab='',xaxt='n',
    #        yaxt='n',xlim=c(0,1),ylim=c(0,1))
    #   legend('left', col = mod_cols[models_plot],lwd=1,
    #          inset=0.03,bty='n',
    #          legend = c('Standard',
    #                     'Non-linear'),
    #          cex=1.1,title = 'Model')
    # } else {
    # draw individual model fit with data
    ind = analysis_data_stan$id==id
    plot(analysis_data_stan$obs_day[ind],
         analysis_data_stan$log_10_vl[ind],
         xlab='', ylab='', 
         xaxt='n', yaxt='n',
         panel.first=grid(), xlim=c(0,7),
         ylim = range(analysis_data_stan$log_10_vl))
    if(counter %% sqrt(K_plots) == 1){
      mtext(text = 'RNA copies per mL',side = 2,
            line = 3,las = 3)
    }
    axis(1, at = c(0,3,7))
    axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                        expression(10^4),
                                        expression(10^6),
                                        expression(10^8)))
    if((counter%%K_plots) >= K_plots - sqrt(K_plots)){
      mtext(text = 'Days',side = 1,line = 2)
    }
    for(mm in models_plot){
      ix = order(analysis_data_stan$obs_day[ind])
      my_xs = analysis_data_stan$obs_day[ind][ix]
      polygon(x = c(my_xs, rev(my_xs)),
              y = c(apply(thetas[[mm]]$preds[,ind],2,
                          quantile,probs=0.025)[ix],
                    rev(apply(thetas[[mm]]$preds[,ind],2,
                              quantile,probs=0.975)[ix])),
              border = NA, 
              col = adjustcolor(mod_cols[mm],alpha.f = .3))
      lines(my_xs,
            colMeans(thetas[[mm]]$preds[,ind])[ix],
            col = mod_cols[mm],lwd=2)
    }
    points(analysis_data_stan$obs_day[ind],
           analysis_data_stan$log_10_vl[ind],pch=16)
    id_map_ind = ID_map$ID_stan==id
    
    mtext(text = paste0(ID_map$ID_key[id_map_ind],
                        '\n',
                        ID_map$Trt[id_map_ind]),
          side = 3, line = -0.5, cex=0.8)
    id=id+1
    # }
    
    counter=counter+1
  }
  
  # add legend to final plot
  plot(NA,NA,xlab='',ylab='',xaxt='n',
       yaxt='n',xlim=c(0,1),ylim=c(0,1))
  legend('left', col = mod_cols[models_plot],lwd=1,
         inset=0.03,bty='n',
         legend = c('Standard',
                    'Non-linear'),
         cex=1.1,title = 'Model')
}


plot_coef_effects = function(stan_out, model_plot, cov_mat, stan_inputs){
  
  thetas = rstan::extract(stan_out[[model_plot]])
  alpha_coefs = apply(thetas$intercept_coefs,2,
                      quantile,probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(alpha_coefs)
  
  cov_names_intercept =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_int[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  cov_names_slope =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_slope[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  plot(alpha_coefs['50%', ], 1:ncol(alpha_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Baseline viral load (fold change)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(alpha_coefs)){
    lines(c(alpha_coefs['10%',i], alpha_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(alpha_coefs['2.5%',i], alpha_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(alpha_coefs), labels = cov_names_intercept,tick = F)
  
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
  
  beta_coefs = apply(thetas$slope_coefs,2, function(x){
    ys = quantile(x, probs=c(0.025,.1,.5,.9,0.975))
    ys = 100*(exp(ys)-1)
    ys
  })
  
  xlims=range(beta_coefs)
  plot(beta_coefs['50%', ], 1:ncol(beta_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Change in rate of clearance (%)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(beta_coefs)){
    lines(c(beta_coefs['10%',i], beta_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(beta_coefs['2.5%',i], beta_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(beta_coefs), labels = cov_names_slope, tick = F)
  x_points = signif(seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = x_points, labels = x_points)
}

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}


checkStrict(make_stan_inputs)
checkStrict(plot_serial_data)
checkStrict(plot_effect_estimates)
checkStrict(plot_individ_data)
checkStrict(make_baseline_table)
checkStrict(plot_coef_effects)
