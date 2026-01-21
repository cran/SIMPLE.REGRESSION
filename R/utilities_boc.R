



# formula check
formula_check <- function(formula, data) {
  
  # check for error in the entered formula
  result_error <- try(as.formula(formula), silent = TRUE)
  if (inherits(result_error, "try-error")) {
    message('\nERROR: The entered formula is not an R formula.')
  } else { formule <- as.formula(formula) }
  
  formula_vars <- all.vars(formule)
  
  # check if all variables are in data
  if (!all(formula_vars %in% colnames(data))) {
    missing_vars <- setdiff(formula_vars, colnames(data))
    message('\nERROR: The following variables were specified in the formula but they do not exist in data.\n')
    print(missing_vars)
  }
  
  DV <- formula_vars[1]
  
  forced <- formula_vars[-1]
  
  output <- list(DV=DV, forced=forced, formule=formule)
  
  return(invisible(output))
}


# formula = 'change ~ music * ticknumb'
# formula = 'change ~ music:ticknumb'
# 
# print(formula_check(formula=formula, data = data_Field_2012))






IV_type_cleanup <- function(df, allIVnoms, verbose=TRUE) {
  
  # converting character IVs to factors
  for (lupe in 1:length(allIVnoms)) {
    if (is.character(df[,allIVnoms[lupe]])) {
      message('\n', allIVnoms[lupe], ' is a character variable. It will be converted into a factor.\n')
      df[,allIVnoms[lupe]] <- as.factor(df[,allIVnoms[lupe]])
    }
  }
  
  # check if all IVs are either numeric or factors
  IVclasses <- sapply(df[,allIVnoms], class)
  if (!all(IVclasses == 'numeric' | IVclasses == 'integer' | IVclasses == 'factor')) {
    message('\n The following variables are not numeric, integers, or factors. Expect errors.')
    names(
      IVclasses[IVclasses != 'numeric' &  IVclasses != 'integer' & IVclasses != 'factor' ])
  }
  
  # test if factor variables have contrasts; if not, add contrasts with names
  IV_factors <- IVclasses[IVclasses == 'factor' ]
  if (length(IV_factors) > 0) {
    for (lupe in 1:length(IV_factors)) {
      
      if (is.null(attr( df[,names(IV_factors[lupe])], "contrasts" ))) {
        
        fac_levels <- levels(df[,names(IV_factors[lupe])])
        
        # print(fac_levels)
        # print(lupe)
        
        # The baseline group is based on alphabetic/numerical order, unless the terms "control"
        # or "Control" or "baseline" or "Baseline" appear in the names of a factor
        # level, in which case that factor level is used as the dummy codes baseline.
        if (length(grep(pattern = 'Control|control|baseline|Baseline', x = fac_levels)) > 0) {
          base_level <- grep(pattern = 'Control|control|baseline|Baseline', x = fac_levels)
        } else {
          base_level <- order(fac_levels, decreasing = FALSE, na.last = TRUE)[1]
        }
        
        custom_contrasts <- contr.treatment(length(fac_levels), base = base_level)
        
        colnames(custom_contrasts) <-
          better_contrast_noms(contrasts(df[,names(IV_factors[lupe])]))
        
        # apply the contrast matrix to the factor variable
        contrasts(df[,names(IV_factors[lupe])]) <- custom_contrasts
      }
    }
  }
  
  # display the contrasts for factors
  if (verbose & length(IV_factors) > 0) {
    for (lupe in 1:length(IV_factors)) {
      message('\nThe contrasts for ', names(IV_factors[lupe]), ' are:\n')
      print(contrasts(df[,names(IV_factors[lupe])]), print.gap=4)
    }
  }
  
  return(invisible(df))
}





descriptives <- function(df, allIVnoms) {
  
  # descriptives for numeric variables
  dfNUM <- df[,allIVnoms]
  dfNUM <- dfNUM[sapply(dfNUM,is.numeric)] # selecting only numeric variables
  if (ncol(dfNUM) != 0) {
    minmax <- t(apply(dfNUM, 2, range))
    descs <- data.frame(Mean=colMeans(dfNUM), SD=apply(dfNUM, 2, sd), 
                        Min=minmax[,1], Max=minmax[,2]) 
    message('\n\nDescriptive statistics for the numeric variables:\n')
    print(round(descs,2), print.gap=4)
  }
  
  # frequencies for factors
  dfFAC <- df[,allIVnoms]
  dfFAC <- dfFAC[sapply(dfFAC,is.factor)] # selecting only factor variables
  if (ncol(dfFAC) != 0) {
    message('\n\nCategory frequencies for the factor variables:\n')
    for (lupe in 1:ncol(dfFAC)) {
      # print(apply((dfFAC), 2, table))
      cat(colnames(dfFAC)[lupe])
      print(table(dfFAC[,lupe]))
    }
  }
}






reg_stats_boc <- function(formule, data, CI_level,
                          prevRsq = NULL, prevModel = NULL, prevmod_lmBF = NULL,
                          MCMC_options) {  # , noms_list) {
  
  model <- lm(formule, data=data, model=TRUE, x=TRUE, y=TRUE)
  
  modelsum <- summary(model)
  
  Rsq <- modelsum$r.squared
  
  # DV name
  DV <- as.character(formule)[[2]] 
  
  # Using the classical R^2 test statistic for (linear) regression designs, 
  # this function computes the corresponding Bayes factor test
  # the Bayes factor (against the intercept-only null)
  Nvars <- ncol(data) - 1
  BF_mod_Rsq <- linearReg.R2stat(R2 = Rsq, N = nrow(data), p = Nvars, 
                                 simple = TRUE)
  
  # Anova Table (Type III tests)
  # only do this if there are no interaction terms in the model
  prednoms <- labels(terms(model))
  contains_Xn <- grepl(":", prednoms)
  anova_table <- NULL
  if (!any(contains_Xn))  anova_table <- ANOVA_TABLE(data=data, model=model) 
  
  # creating a new version of the raw data, from the lm output, so that can 
  # provide stats for dummy codes
  # modeldata is needed for moderation analyses (next)
  modeldata <- data.frame(model$y, model$x[,2:ncol(model$x)])
  colnames(modeldata) <- c(DV,colnames(model$x)[2:ncol(model$x)])
  
  modelsum$coefficients <- cbind(modelsum$coefficients, confint(model))
  
  partialRcoefs <- 
    PARTIAL_COR(data=cor(modeldata), Y=DV, X=names(model$coefficients)[-1], 
                Ncases=nrow(data), verbose=FALSE)
  
  partialRcoefs <- cbind(partialRcoefs$betas, partialRcoefs$Rx_y, 
                         partialRcoefs$R_partials, partialRcoefs$R_semipartials)
  colnames(partialRcoefs) <- c('beta','r','partial.r','semipartial.r')
  
  
  modeldata$predicted <- model$fitted.values
  
  # casewise diagnostics
  modeldata$residuals <- resid(model)
  modeldata$residuals_standardized <- rstandard(model)
  modeldata$residuals_studentized <- rstudent(model)
  modeldata$cooks_distance <- cooks.distance(model)
  modeldata$dfbeta <- dfbeta(model)
  modeldata$dffit <- dffits(model)
  modeldata$leverage <- hatvalues(model)
  modeldata$covariance_ratios <- covratio(model)
  
  collin_diags <- Collinearity(model.matrix(model), verbose=FALSE)
  
  # Bayes factors - for the predictors, from just t & df
  Bayes_BFs <- BF_reg_preds(mod_sum = modelsum)
  
  # Bayes HDIs
  Bayes_HDIs <- chains <- autocorrels <- SDs <- NULL
  
  if (MCMC_options$MCMC)  {
    mod_HDIs <- Bayes_HDIs_reg_preds(formule = formule, 
                                     data = data, 
                                     CI_level = CI_level,
                                     MCMC_options = MCMC_options) #,  noms_list = noms_list)
    
    Bayes_HDIs  <- mod_HDIs$Bayes_HDIs
    
    chains      <- mod_HDIs$chains
    
    autocorrels <- mod_HDIs$autocorrels
    
    # mod_lmBF    <- mod_HDIs$mod_lmBF
    # SDs         <- mod_HDIs$SDs
  }
  
  mod_lmBF <- lmBF(formule, data, columnFilter="ID") 
  
  # current vs previous model comparisons	
  Rsqch <- fsquared <- fish <- BF_change <- NULL
  if (!is.null(prevRsq) & !is.null(prevModel) & !is.null(prevmod_lmBF)) {
    
    # lm model comparisons
    Rsqch <- Rsq - prevRsq
    fsquared <- (Rsq - prevRsq) / (1 - Rsq)	
    fish <- anova(model, prevModel)
    
    # lmBF model comparisons
    prevmod_lmBF@data <- mod_lmBF@data
    BF_change <- mod_lmBF / prevmod_lmBF
    BF_change <- as.numeric(unlist(extractBF(BF_change))[1])
    # BF_interps(BF_M2 = BF_change)
  }
  
  # # current vs previous lmBF model comparisons	
  # BF_change <- NULL
  # if (MCMC_options$MCMC & !is.null(prevmod_lmBF)) {
  # 
  #   prevmod_lmBF@data <- mod_lmBF@data
  # 
  #   BF_change <- mod_lmBF / prevmod_lmBF
  # 
  #   BF_change <- as.numeric(unlist(extractBF(BF_change))[1])
  # 
  #   BF_interps(BF_M2 = BF_change)
  # }
  
  cormat <- cor(modeldata[,c(DV,names(model$coefficients)[-1])])
  
  output <- list(model=model, modelsum=modelsum, Rsq=Rsq, anova_table=anova_table, 
                 Rsqch=Rsqch, fsquared=fsquared, fish=fish,
                 modeldata=modeldata, partialRcoefs=partialRcoefs, 
                 collin_diags=collin_diags, BF_mod_Rsq=BF_mod_Rsq,
                 Bayes_BFs=Bayes_BFs, Bayes_HDIs=Bayes_HDIs, 
                 mod_lmBF=mod_lmBF, BF_change=BF_change,
                 chains=chains, autocorrels=autocorrels,  # , SDs=SDs
                 MCMC_options=MCMC_options, cormat=cormat) # , noms_list=noms_list)
  
  return(invisible(output))
}







# Bayes factors - for the predictors, from just t & df from an lm model summary
BF_reg_preds <- function(mod_sum) {
  
  t_values <- mod_sum$coefficients[,3]
  BF_preds_h10 <- unlist(sapply(t_values, ttest.tstat, 
                                n1 = (mod_sum$df[2] + 1), 
                                simple=TRUE))
  # using df.residual(model) + 1 for the n1 value in ttest.tstat bec
  # assuming that the function computes the df by n1 - 1
  BF_preds_h01 <- 1 / BF_preds_h10
  Bayes_BFs <- data.frame(cbind(mod_sum$coefficients[,1], t_values, BF_preds_h10, BF_preds_h01))
  colnames(Bayes_BFs) <- cbind('b', 't', 'BF_10', 'BF_01')
  Bayes_BFs$Interpretation <- apply( cbind(Bayes_BFs$BF_10, Bayes_BFs$BF_01), 1, BF_interps_2)
  
  return(invisible(Bayes_BFs))
}




# Bayes HDIs
Bayes_HDIs_reg_preds <- function(formule, data, CI_level = 95, 
                                 MCMC_options, # , noms_list, 
                                 model_type = 'OLS', famille = NULL) {
  
  # 2 functions/packages are used for the HDIs
  
  # BayesFactor package -- PROB = factor contrast names, but it does provide
  # for model comparisons (via model division, & not by a function)
  
  # rstanarm::stan_lm provides the same contrasts, with the same names, for
  # factors as the lm function, but no way of comparing models
  # the loo package is used for these model comparisons, but it has many
  # dependencies & runs leave-one-out analyses
  
  # so below, when there are no factors, only the BayesFactor package is used
  # when there are factors in a model, rstanarm::stan_lm is used for the
  # coefs & chains, & the BayesFactor package is also used for the model
  # comparisons. This means that the MCMC analyses are run twice (yuk)
  
  
  # # BayesFactor package -- PROB = factor contrast names, +
  # mod_lmBF <- lmBF(formule, data, columnFilter="ID") 
  # 
  # chains <- posterior(mod_lmBF, iterations = MCMC_options$Nsamples, 
  #                     thin = MCMC_options$thin, progress = FALSE)
  
  # rstanarm
  # mod_stan <- rstanarm::stan_lm(formula = formule, data = data,
  #                                  prior = R2(0.5), seed = 12345,
  #                                  chains = 1, iter = MCMC_options$Nsamples)
  
  if (model_type == 'OLS') {
    mod_stan <- suppressWarnings(
      rstanarm::stan_lm(formula = formule, data = data,
                                  prior = NULL, seed = 12345,
                                  chains = 4, iter = MCMC_options$Nsamples,
                                  thin = MCMC_options$thin,
                                  warmup = MCMC_options$burnin,
                                  refresh = 0)) #, control = list(adapt_delta = 0.99)) )
  }
  
  if (model_type == 'LOGISTIC') {
    
    # rstanarm: The error "All outcome values must be 0 or 1 for Bernoulli models" 
    # occurs because the rstanarm package requires the response variable 
    # in a Bernoulli model to be numeric and strictly coded as 0 for 
    # failure and 1 for success. The model cannot handle other numeric 
    # values (like proportions or counts greater than 1) or factor/character 
    # data types in the outcome variable. 
    
    # convert DV to numeric, 0 or 1
    dd <- ifelse( data[,1] == levels(data[,1])[1], 0, 1)
    data[,1] <- dd        
    # is.numeric(data[,1]); print(table(data[,1]))
    
    mod_stan <- rstanarm::stan_glm(formule, data = data, 
                                   family = "binomial",
                                   refresh = 0, algorithm="sampling", 
                                   chains = 4, iter = MCMC_options$Nsamples,
                                   thin = MCMC_options$thin,
                                   warmup = MCMC_options$burnin)
    # print(prior_summary(mod_stan))
  }
  
  if (model_type == 'COUNT') {
    
    mod_stan <- rstanarm::stan_glm(formule, data = data, 
                                   family = famille,
                                   refresh = 0, algorithm="sampling", 
                                   chains = 4, iter = MCMC_options$Nsamples,
                                   thin = MCMC_options$thin,
                                   warmup = MCMC_options$burnin)
    # print(prior_summary(mod_stan))
  }
  
  MCMC_mod_coefs <- cbind( coef(mod_stan), 
                           rstanarm::posterior_interval(mod_stan, prob=CI_level * .01)
                           [1:length(coef(mod_stan)),])
  # round(MCMC_mod_coefs, 3)
  
  # extract samples as a data frame
  chains <- as.data.frame(mod_stan) # ;  dim(chains_df); head(chains_df)
  
  # # remove burn-in period
  # chains <- chains[(MCMC_options$burnin+1):nrow(chains),]
  
  colnames(chains)[colnames(chains) == "(Intercept)"] <- "Intercept"
  
  # removing unneeded columns from OLS chains
  if (model_type == 'OLS') {
    colnames(chains)[colnames(chains) == "log-fit_ratio"] <- "log_fit_ratio"
    chains$sigma <- NULL
    chains$log_fit_ratio <- NULL
    chains$R2 <- NULL
  }
  
  Nests <- ncol(chains)
  
  # autocorrelations
  lag.max <- 10
  autocorrels <- matrix(NA, (lag.max+1), ncol(chains))
  for (lupe in 1:ncol(chains)) 
    autocorrels[,lupe] <- as.matrix(acf(chains[,lupe], type = "correlation", 
                                        plot = FALSE, lag.max = lag.max)$acf)
  colnames(autocorrels) <- colnames(chains)
  rownames(autocorrels) <- paste0("Lag ", 0:(nrow(autocorrels)-1))
  # colnames(autocorrels)[colnames(autocorrels) == "(Intercept)"] <- "Intercept"
  # colnames(autocorrels)[colnames(autocorrels) == "log-fit_ratio"] <- "log_fit_ratio"
  # autocorrels <- subset(autocorrels, select = -c(sigma, log_fit_ratio, R2))
  round(autocorrels,4)
  
  # removing the COVARS & other less important columns
  # keypreds <- c(noms_list$IV, noms_list$MOD_terms, noms_list$PROD)
  # keepers <- intersect(colnames(autocorrels), keypreds)
  # autocorrels <- autocorrels[,keepers]
  
  # parameters estimates & CIs
  ests <- colMeans(chains[,1:Nests])
  quant_size <- (100 - CI_level) / 2
  quant_lb <- quant_size * .01
  quant_ub <- (100 - quant_size) * .01
  ests_ci_lb <- apply(chains[,1:Nests],2,quantile,probs=quant_lb)
  ests_ci_ub <- apply(chains[,1:Nests],2,quantile,probs=quant_ub)
  Bayes_HDIs <- cbind(ests, ests_ci_lb, ests_ci_ub)
  Bayes_HDIs <- Bayes_HDIs[rownames(Bayes_HDIs) != "sigma", ]
  
  
  if (model_type == 'OLS') {
    
    # standardized estimates & CIs  -- OLD
    # first get the names of the variables in formule
    # varnoms <- all.vars(formule)
    # SDs <- apply(data[,varnoms],2,sd)
    # ests_z <- c(NA, (ests[-1] * SDs[-1]) / SDs[1])
    # ests_z_ci_lb <- c(NA, (ests_ci_lb[-1] * SDs[-1]) / SDs[1])
    # ests_z_ci_ub <- c(NA, (ests_ci_ub[-1] * SDs[-1]) / SDs[1])
    
    # standardized estimates & CIs
    # if a variable is in the list of continuous predictors, compute beta using
    # the usual formula
    # if a variable is NOT in the list of continuous predictors, compute beta using
    # only the SD of the DV (for partial standardization)
    
    varnoms <- all.vars(formule)[-1]
    IV_types <- sapply(data[,varnoms], class)
    contin_vars <- IV_types[IV_types == 'numeric' | IV_types == 'integer']
    ests_z <- matrix(NA, nrow(Bayes_HDIs), 3)
    for (lupe in 2:nrow(Bayes_HDIs)) {
      
      if (rownames(Bayes_HDIs)[lupe] %in% names(contin_vars)) {
        ests_z[lupe,] <- (Bayes_HDIs[lupe,] * sd(data[,rownames(Bayes_HDIs)[lupe]])) / sd(data[,1])
        
      } else {
        ests_z[lupe,] <- (Bayes_HDIs[lupe,] * 1) / sd(data[,1])
      }
    }
    Bayes_HDIs <- cbind(Bayes_HDIs, ests_z)
    rownames(Bayes_HDIs)[1] <- 'Intercept'
    colnames(Bayes_HDIs) <- cbind('b', 'b_ci_lb', 'b_ci_ub',
                                  'beta', 'beta_ci_lb', 'beta_ci_ub')
    round(Bayes_HDIs,5)
  }
  
  if (model_type == 'LOGISTIC' | model_type == 'COUNT') {
    
    # providing the raw coefs & the corresponding exponentiated coefs (not std coefs)
    Bayes_HDIs <- cbind(Bayes_HDIs, exp(Bayes_HDIs))
    
    if (model_type == 'LOGISTIC')
      colnames(Bayes_HDIs) <- c('b', 'b_ci_lb', 'b_ci_ub', 'OR','OR_ci_lb','OR_ci_ub')
    
    if (model_type == 'COUNT')
      colnames(Bayes_HDIs) <- c('b', 'b_ci_lb', 'b_ci_ub', 
                                'IRR', 'IRR_ci_lb', 'IRR_ci_ub')  
  }
  
  
  # https://cran.r-project.org/web/packages/BayesFactor/vignettes/manual.html#regression
  # The results are quite similar (to OLS coefs), apart from the intercept. This is due to the 
  # Bayesian model centering the covariates before analysis, so the mu parameter 
  # is the mean of  rather than the expected value of the response variable 
  # when all uncentered covariates are equal to 0.
  
  output <- list(Bayes_HDIs=Bayes_HDIs, chains=chains, 
                 autocorrels=autocorrels)  # SDs=SDs, 
  
  return(invisible(output))
}   





Bayes_HDI_plot <- function(chain_dat, Bayes_HDIs, CI_level, 
                           HDI_plot_est_type = 'standardized') {   # , SDs
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(ask=TRUE)
  
  # if (HDI_plot_est_type == 'standardized') {
  #   # chain_dat_z <- matrix(NA, nrow(chain_dat), ncol(chain_dat))
  #   for (lupe in 1:ncol(chain_dat))  chain_dat[,lupe] = (chain_dat[,lupe] * SDs[lupe+1]) / SDs[1]
  # }
  colMeans(chain_dat)
  
  for (lupe in 1:ncol(chain_dat)) {
    
    if (HDI_plot_est_type == 'raw') {
      HDIint <- Bayes_HDIs[lupe, c(2:3)]
      
      xlab = paste("raw b for ", colnames(chain_dat)[lupe], sep='')
    }
    
    # if (HDI_plot_est_type == 'standardized') {
    #   HDIint <- Bayes_HDIs[lupe, c(5:6)]
    # 
    #   xlab = paste("Beta for ", colnames(chain_dat)[lupe], sep='')
    # }
    
    dendat = density(chain_dat[,lupe])
    
    plot(density(chain_dat[,lupe]),
         main=paste("Bayesian Posterior Distribution for\n", colnames(chain_dat)[lupe]), 
         font.main = 1.5, cex.main=1.5,
         xlab = xlab,
         col = "skyblue",
         # ylim = c(0,ceiling(max(dendat$y))),
         ylim = c(0,  (max(dendat$y))),
         xlim = c(min(chain_dat[,lupe]), max(chain_dat[,lupe])), 
         cex.lab = 1.2)
    
    polygon(cbind(c(0, dendat$x, max(dendat$x)), c(0, dendat$y, min(dendat$y))), col="skyblue")
    
    abline(v=0, lty="dotted", lwd=1, col="red")
    abline(v=HDIint[1], lty="dotted", lwd=1.5, col="black")
    abline(v=HDIint[2], lty="dotted", lwd=1.5, col="black")
    
    lines(HDIint, c(0,0), lwd=10, lend=1)
    
    text( HDIint[1], 0, bquote(.(signif(HDIint[1],3))), adj=c((0.7+.4),-15.5), cex=1.0)
    text( HDIint[2], 0, bquote(.(signif(HDIint[2],3))), adj=c((0.7-.9),-15.5), cex=1.0)
    
    text( mean(HDIint), 0, bquote(.(CI_level) * "% HDI"), adj=c(.5,-0.8), cex=1.0)
  }
} 





BF_interps <- function(BF_10 = NULL, BF_01 = NULL, BF_M1 = NULL, BF_M2 = NULL) {
  
  # null & alternative
  if (is.null(BF_10) & !is.null(BF_01))  BF_10 <- 1 / BF_01
  if (is.null(BF_01) & !is.null(BF_10))  BF_01 <- 1 / BF_10
  
  if (!is.null(BF_10) & !is.null(BF_01)) {
    
    message('\nBayes Factor for the null vs. alternative model: ', round(BF_01,2))
    message('\nBayes Factor for the alternative vs. null model: ', round(BF_10,2))
    
    message('\nThe Bayes factor suggests that the observed data are ', round(BF_01,2), ' times more likely')
    message('to occur under the null hypothesis than under the alternative hypothesis.')
    message('\nAlternatively, the observed data are ', round(BF_10,2), ' times more likely to occur')
    message('under the alternative hypothesis than under the null hypothesis.')
    
    # 2016 - Diagnosing the Misuse of the Bayes Factor in Applied Research p 4
    if (BF_01 > BF_10) {
      message('\nSomeone with no prior preference for either hypothesis (i.e., prior odds = 1)')
      message('should now believe that the null model is ', round(BF_01,2), ' times more probable')
      message('than the alternative model (i.e., posterior odds = ', round(BF_01,2), ' x 1 = ', round(BF_01,2), ').')
    }
    
    if (BF_10 > BF_01) {
      message('\nSomeone with no prior preference for either hypothesis (i.e., prior odds = 1)')
      message('should now believe that the alternative model is ', round(BF_10,2), ' times more probable')
      message('than the null model (i.e., posterior odds = ', round(BF_10,2), ' x 1 = ', round(BF_10,2), ').\n')
    }
    
    if (BF_10 > BF_01) hyp <- 'alternative'
    if (BF_01 > BF_10) hyp <- 'null'
    
    # Jeffreys' scale | Grades or categories of evidence for the Bayes factor
    # Lee M. D. and Wagenmakers, E.J. (2014) Bayesian cognitive modeling: A 
    # practical course, Cambridge University Press.
    maxBF <- max(c(BF_10, BF_01))
    if (maxBF >= 1  & maxBF < 3)   descr <- 'Anecdotal'
    if (maxBF >= 3  & maxBF < 10)  descr <- 'Moderate'
    if (maxBF >= 10 & maxBF < 30)  descr <- 'Strong'
    if (maxBF >= 30 & maxBF < 100) descr <- 'Very Strong'
    if (maxBF >  100)              descr <- 'Extreme'
    
    message('This level of evidence in favor of the ', hyp, ' hypothesis is') 
    message('considered "', descr,'" according to Lee & Wagenmakers (2014).')
  }
  
  # M1 & M2
  if (is.null(BF_M2) & !is.null(BF_M1))  BF_M2 <- 1 / BF_M1
  if (is.null(BF_M1) & !is.null(BF_M2))  BF_M1 <- 1 / BF_M2
  
  if (!is.null(BF_M2) & !is.null(BF_M1)) {
    
    message('\n    Bayes Factor for the previous vs. current models: ', round(BF_M1,2))
    message('\n    Bayes Factor for the current vs. previous models: ', round(BF_M2,2))
    
    message('\n    The Bayes factor suggests that the data are ', round(BF_M1,2), ' times more likely')
    message('    to occur under the previous model than under the current model.')
    message('\n    Alternatively, the data are ', round(BF_M2,2), ' times more likely to occur')
    message('    under the current model than under the previous model.')
    
    if (BF_M2 > BF_M1) hyp <- 'current'
    if (BF_M1 > BF_M2) hyp <- 'previous'
    
    maxBF <- max(c(BF_M2, BF_M1))
    if (maxBF >= 1  & maxBF < 3)   descr <- 'Anecdotal'
    if (maxBF >= 3  & maxBF < 10)  descr <- 'Moderate'
    if (maxBF >= 10 & maxBF < 30)  descr <- 'Strong'
    if (maxBF >= 30 & maxBF < 100) descr <- 'Very Strong'
    if (maxBF >  100)              descr <- 'Extreme'
    
    message('\n    This level of evidence in favor of the ', hyp, ' hypothesis is') 
    message('    considered "', descr,'" according to Lee & Wagenmakers (2014).')
  }
}

# BF_interps(BF_10 = 3.2, BF_01 = 1/3.2)
# 
# BF_interps(BF_10 = 3.2)
# 
# BF_interps(BF_M1 = 3.2, BF_M2 = 1/3.2)




# BF_interps_2 <- function(BF_M1, BF_M2) {  
BF_interps_2 <- function(dat) {
  
  BF_M1 <- dat[1]
  BF_M2 <- dat[2]
  
  if (BF_M1 > BF_M2) hyp <- 'H1'
  if (BF_M2 > BF_M1) hyp <- 'H0'
  
  # Jeffreys' scale | Grades or categories of evidence for the Bayes factor
  # Lee M. D. and Wagenmakers, E.J. (2014) Bayesian cognitive modeling: A 
  # practical course, Cambridge University Press.
  maxBF <- max(c(BF_M1, BF_M2))
  if (maxBF >= 1  & maxBF < 3)   descr <- paste('Anecdotal', hyp, sep = ' evidence for ')
  if (maxBF >= 3  & maxBF < 10)  descr <- paste('Moderate', hyp, sep = ' evidence for ')
  if (maxBF >= 10 & maxBF < 30)  descr <- paste('Strong', hyp, sep = ' evidence for ')
  if (maxBF >= 30 & maxBF < 100) descr <- paste('Very Strong', hyp, sep = ' evidence for ')
  if (maxBF >  100)              descr <- paste('Extreme', hyp, sep = ' evidence for ')
  
  return(invisible(descr))
}





VIFtol_messages <- function(VIFtol) {
  message('\nThe mean Variance Inflation Factor (VIF) = ', round(mean(VIFtol),3))
  message('\nMulticollinearity is said to exist when VIF values are > 5, or when Tolerance values')
  message('are < 0.1. Multicollinearity may be biasing a regression model, and there is reason')
  message('for serious concern, when VIF values are greater than 10.\n\n')
}

CondInd_messages <- function(CondInd) {
  message('\n')
  print(round(CondInd,3), print.gap=4)
  message('\nThe coefficients for the Intercept and for the predictors are the Variance Proportions.')
  message('\nEigenvalues that differ greatly in size suggest multicollinearity. Multicollinearity is')
  message('said to exist when the largest condition index is between 10 and 30, and evidence for')
  message('the problem is strong when the largest condition index is > 30. Multicollinearity is')
  message('also said to exist when the Variance Proportions (in the same row) for two variables are')
  message('above .80 and when the corresponding condition index for a dimension is higher than 10 to 30.\n\n')
}





# display output from reg_stats_boc
resultats <- function(xx) {
  
  message('\n\nmultiple R = ', round(sqrt(xx$Rsq),3), 
          '   multiple R-squared = ', round(xx$Rsq,3),
          '   adjusted R-squared = ', round(xx$modelsum$adj.r.squared,3))
  
  Fstats <- xx$modelsum$fstatistic
  pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
  
  message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2],  
          '   df_denom = ', Fstats[3], '   p-value = ', round(pvalue,6), '\n')
  
  BF_interps(BF_10 = xx$BF_mod_Rsq, BF_01 = (1 / xx$BF_mod_Rsq) )
  
  
  if (!is.null(xx$Rsqch)) {
    # current vs previous lm model comparisons	
    message('\n\nCurrent vs previous model comparison:')
    message('\n    Rsquared change = ', round(xx$Rsqch,3))
    message('\n    f-squared change = ', round(xx$fsquared,3))
    
    message('\n    F = ', round(xx$fish$F[2],2), '   df_num = ', 1,  
            '   df_denom = ', xx$fish$Res.Df[1], '   p-value = ', 
            round(xx$fish$'Pr(>F)'[2],6), '\n')
    
    # current vs previous lmBF model comparisons	
    if (!is.null(xx$BF_change)) BF_interps(BF_M2 = xx$BF_change)
  }
  
  if (!is.null(xx$anova_table)) {
    message('\n\nAnova Table (Type III tests):\n')
    print(round_boc(xx$anova_table,3), print.gap=4)
  }
  
  message('\nModel coefficients:\n')
  print(round_boc(xx$modelsum$coefficients,3), print.gap=4)
  
  message('\n\nBayes Factors for the predictors:\n')
  print(round_boc(xx$Bayes_BFs, round_non_p = 2, round_p = 5), print.gap=4)
  
  message('\n\nBeta, r, partial correlations, & semi-partial correlations:\n')
  print(round(xx$partialRcoefs,3), print.gap=4)	
  
  if (xx$MCMC_options$MCMC) {
    
    message('\n\nBayesian MCMC chains:')
    message('\n    Number iterations (samples): ', xx$MCMC_options$Nsamples)
    message('\n    Thinning interval:  ', xx$MCMC_options$thin)
    message('\n    Burn-in period:  ', xx$MCMC_options$burnin)
    message('\n    Final chain length (# of samples): ', nrow(xx$chains))
    # message('\n    The priors were the default BayesFactor::lmBF function priors')
    message('\n    The priors were the default rstanarm::stan_lm function priors')
    
    message('\n\nAutocorrelations for the MCMC chains:\n')
    print(round(xx$autocorrels,3), print.gap=4)
    
    message('\n\nBayesian HDIs for the raw and standardized coefficients:\n')
    print(round_boc(xx$Bayes_HDIs, round_non_p = 2, round_p = 5), print.gap=4)
    
    # notice of partial standardization of factor contrast effects
    # were there any predictors in the model factors?
    factor_vars <- names(xx$model$xlevels)
    if (length(factor_vars) > 0) {
      
      # vector of term names that were not for the continuous variables
      pred_vars <- names(attr(xx$model$terms, "dataClasses"))[-1]
      contin_vars <- setdiff(pred_vars, factor_vars)
      all_terms <- names(xx$modeldata)[-1]
      # names of the regression diagnostic variables -- for exclusion from all_terms
      nomsdiags <- c("residuals_standardized", "residuals_studentized", "cooks_distance",
                     "dfbeta", "dffit", "leverage", "covariance_ratios", "predicted", "residuals")
      non_contin_terms <- setdiff(all_terms, c(contin_vars, factor_vars , nomsdiags))
      
      message('\nNote. The above beta coefficients for the terms listed below were only')
      message('      partially standardized. Specifically, only the DV was placed in')
      message('      standard deviation metric, to facilitate interpretations:')
      
      cat('\t', non_contin_terms, sep='\n\t')
    }
  }
  
  message('\n\nCorrelation matrix:\n')
  print(round(xx$cormat,3), print.gap=4)		
  
  message('\n\nCollinearity Diagnostics:')
  VIFtol_messages(xx$collin_diags$VIFtol[,2])
  CondInd_messages(xx$collin_diags$CondInd)
}





# better names for dummy contrasts function
better_contrast_noms <- function(contrs, nom_max_len = 21) {
  
  new_noms <- c()
  
  # the baseline (all 0s)
  base_level <- apply(contrs, 1, function(row) all(row == 0)) 
  base_level <- rownames(contrs)[which(base_level == TRUE)]
  base_level <- gsub(" ", "_", base_level)
  base_level <- substr(base_level, 1, 10)
  if (substr(base_level, nom_max_len, nom_max_len) == '_')  
    base_level <- substr(base_level, 1, (nom_max_len - 1))
  
  new_nom_max_len <- nom_max_len - nchar(base_level) - 4
  
  for (lupe in 1:ncol(contrs)) {
    
    nom_max_len <- 20
    
    comp <- names(which(contrs[,lupe] == 1))
    comp <- gsub(" ", "_", comp)
    comp <- substr(comp, 1, new_nom_max_len)
    if (substr(comp, new_nom_max_len, new_nom_max_len) == '_')  
      comp <- substr(comp, 1, (new_nom_max_len - 1))
    
    # new_noms <- c(new_noms, paste('__', paste(comp, base_level, sep='_vs_'), sep=''))
    new_noms <- c(new_noms, paste(comp, base_level, sep='_vs_'))
  }
  return(invisible(new_noms))
}




# finds an empty space for locating the plot legend
legend_loc <- function(x, y, xlim, ylim) {
  
  # x = c(0,1);  xf = cut(x, breaks=3)
  # x = seq(0, 1, by=.05)
  # xlim = c(0,1)
  
  x_thirds <- (abs(xlim[1] - xlim[2])) / 3
  
  x_breaks <- c( xlim[1], (xlim[1] + x_thirds), (xlim[1] + (x_thirds * 2)), xlim[2] )
  
  xf <- cut(x, breaks = x_breaks, include.lowest = TRUE)
  
  y_thirds <- (abs(ylim[1] - ylim[2])) / 3
  
  y_breaks <- c( ylim[1], (ylim[1] + y_thirds), (ylim[1] + (y_thirds * 2)), ylim[2] )
  
  yf <- cut(y, breaks = y_breaks, include.lowest = TRUE)
  
  tab <- table(yf, xf);  #print(tab)
  
  # rearrange to match the R positions
  tab <- rbind(tab[3,], tab[2,], tab[1,])
  
  tabmin <- min(tab)
  
  if        (tab[1,3] == tabmin) { loc <- 'topright'
  } else if (tab[1,1] == tabmin) { loc <- 'topleft'
  } else if (tab[1,2] == tabmin) { loc <- 'top'
  } else if (tab[3,3] == tabmin) { loc <- 'bottomright'
  } else if (tab[3,1] == tabmin) { loc <- 'bottomleft'
  } else if (tab[3,2] == tabmin) { loc <- 'bottom'
  } else if (tab[2,3] == tabmin) { loc <- 'right'
  } else if (tab[2,1] == tabmin) { loc <- 'left'
  } else if (tab[2,2] == tabmin) { loc <- 'center'
  }
  
  return(invisible(loc))
}





plotfun <- function(testdata, list_xlevels, DV_predicted, CIs, kind,
                    xlim, ylim, xlab, ylab, title, cols_user,
                    IV_focal_1, IV_focal_1_values, IV_focal_2, IV_focal_2_values) {
  
  # if IV_focal_1 is NOT a factor
  if (!IV_focal_1 %in% names(list_xlevels)) {
    
    if (is.null(xlim))  xlim <- c(min(testdata[,IV_focal_1]), max(testdata[,IV_focal_1]))
    
    if (is.null(IV_focal_2)) {
      
      plot(testdata[,DV_predicted] ~ testdata[,IV_focal_1], type = 'n', 
           xlim = xlim, 
           xlab = xlab,
           ylim = ylim, 
           ylab = ylab,
           main = title)
      
      if (CIs) {
        polygon(c(rev(testdata[,IV_focal_1]), testdata[,IV_focal_1]), 
                c(rev(testdata$ci_ub), testdata$ci_lb), col = 'grey95', border = NA)
        lines(testdata[,IV_focal_1], testdata$ci_ub, col = 'red', lwd = .4)  # lty = 'dashed',
        lines(testdata[,IV_focal_1], testdata$ci_lb, col = 'red', lwd = .4)
      }
      lines(testdata[,IV_focal_1], testdata[,DV_predicted],  col = 'black', lwd = 1.3)
      
      grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
    }
    
    if (!is.null(IV_focal_2)) { 
      
      # colours for the IV_focal_2_values
      if (length(cols_user) >= length(IV_focal_2_values))  {
        cols <- cols_user[1:length(IV_focal_2_values)]
      } else {cols <- rainbow(length(IV_focal_2_values))}
      
      dum <- testdata[testdata[IV_focal_2] == IV_focal_2_values[1],]
      
      plot(dum[,DV_predicted] ~ dum[,IV_focal_1], type = 'l', 
           xlim = xlim, 
           xlab = xlab,
           ylim = ylim, 
           ylab = ylab,
           main = title,
           col = cols[1])
      
      for (lupe in 2:length(IV_focal_2_values)) {
        
        dum <- testdata[testdata[IV_focal_2] == IV_focal_2_values[lupe],]
        
        lines(dum[,IV_focal_1], dum[,DV_predicted], col = cols[lupe], lwd = 1.3)
      }
      
      legvalues <- IV_focal_2_values
      if (is.numeric(legvalues))  legvalues <- round(legvalues,3)
      
      leg_location <- legend_loc(x = testdata[,IV_focal_1], y = testdata[,DV_predicted], 
                                 xlim = xlim, ylim = ylim)
      
      legend(leg_location, legend=legvalues, title=IV_focal_2, col=cols,  
             bty="n", lty=1, lwd=2, cex = 0.75)   
      
      grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
    }
  } 
  
  
  # if IV_focal_1 is a factor
  if (IV_focal_1 %in% names(list_xlevels)) {
    
    # # colours for the IV_focal_1_values
    # if (length(cols_user) >= length(IV_focal_1_values))  {
    #   cols <- cols_user[1:length(IV_focal_1_values)]
    # } else {cols <- rainbow(length(IV_focal_1_values))}
    
    # colours for the IV_focal_2_values
    if (length(cols_user) >= length(IV_focal_2_values))  {
      cols <- cols_user[1:length(IV_focal_2_values)]
    } else {cols <- rainbow(length(IV_focal_2_values))}
    
    # if (is.null(ylim)) {
    #   if (CIs) {ylim <- range( pretty(testdata$ci_lb), pretty(testdata$ci_ub))
    #   } else {ylim <- range( pretty(testdata[,DV_predicted]))}
    # }
    # 
    # if (kind != 'LOGISTIC')  ylim[2] <- ylim[2] + (ylim[2] * .05)
    
    ylim <- better_ylim(ylim, buffer = 0.1)
    
    # to avoid an upside down bar plot when all y values are negative, 
    # reverse the Y-axis direction (values go from high to low)
    if (all(ylim < 0))  ylim <- rev(range(ylim))
    
    # functions to produce wrapped group names for the X-axis
    # https://stackoverflow.com/questions/20241065/r-barplot-wrapping-long-text-labels
    # Core wrapping function
    wrap.it <- function(x, len)
    { sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE) }
    # Call this function with a list or vector
    wrap.labels <- function(x, len) {
      if (is.list(x)) { lapply(x, wrap.it, len)
      } else { wrap.it(x, len) }
    }
    wr.lap <- wrap.labels(IV_focal_1_values, 8)
    
    
    if (is.null(IV_focal_2)) {  
      
      don <- t(testdata[DV_predicted]); colnames(don) <- IV_focal_1_values
      
      plot_bar <- barplot(don, beside=TRUE, legend.text=FALSE, col=cols,
                          yaxp = c(min(ylim), max(ylim), 4),
                          ylim= rev(range(ylim)), #ylim, #yaxt="n", 
                          xpd=FALSE,
                          xlab = IV_focal_1,
                          ylab = ylab,
                          main = title,
                          names.arg = wr.lap)
      
      if (CIs) {
        arrows(plot_bar, t(testdata$ci_ub), plot_bar, t(testdata$ci_lb), 
               lwd = 1.0, angle = 90, code = 3, length = 0.05)
      }
      
      # graphics::legend("top", legend = IV_focal_1_values, bty="n", lwd=4, col=cols, cex = .80)
    }
    
    if (!is.null(IV_focal_2)) {  
      
      plotdon_DV  <- unlist(testdata[DV_predicted])
      
      plotdon_IV1 <- unlist(testdata[IV_focal_1])
      
      plotdon_IV2 <- unlist(testdata[IV_focal_2])
      
      # if IV_focal_2 is not a factor, round to make the ticks shorter
      if (!IV_focal_2 %in% names(list_xlevels))  plotdon_IV2 <- round(plotdon_IV2,2)
      
      plot_bar <- barplot(plotdon_DV ~ plotdon_IV2 + plotdon_IV1, beside=TRUE,
                          legend.text=FALSE, col=rep(cols, length(IV_focal_2_values)), 
                          yaxp = c(min(ylim), max(ylim), 4),
                          ylim = ylim, 
                          xpd=FALSE,
                          xlab = IV_focal_1,
                          ylab = ylab,
                          main = title,
                          names.arg = wr.lap)
      
      if (CIs) {
        arrows(plot_bar, t(testdata$ci_ub), plot_bar, t(testdata$ci_lb), 
               lwd = 1.0, angle = 90, code = 3, length = 0.05)
      }
      graphics::legend("top", legend = IV_focal_2_values, title = IV_focal_2, 
                       bty="n", lwd=4, col=cols, cex = .80)
    }
  }

}





# function to expand (y-axis plot) limits by a percentage (e.g., 10%)
better_ylim <- function(y, buffer = 0.1) {
  r <- range(y, na.rm = TRUE)
  diff <- diff(r)
  return(c(r[1] - diff * buffer, r[2] + diff * buffer))
}






produce_ylim <- function(data_vector, buffer_percent = 0.1) {
  # Remove NA and infinite values
  clean_data <- data_vector[is.finite(data_vector)]
  
  if (length(clean_data) == 0) {
    warning("No finite data points to determine range.")
    return(NULL)
  }
  
  data_range <- range(clean_data)
  buffer_val <- buffer_percent * diff(data_range)
  
  lower_bound <- data_range[1] - buffer_val
  upper_bound <- data_range[2] + buffer_val
  
  return(c(lower_bound, upper_bound))
}

# Example usage:
my_data <- c(1, 5, 10, 15, 20)
ideal_limits <- produce_ylim(my_data, buffer_percent = 0.05)
plot(my_data, ylim = ideal_limits)





# Goodness of Fit
GoF_stats <- function(model) {
  
  # X2_Pearson <- sum(residuals(model, type = "pearson")^2)
  
  # https://stackoverflow.com/questions/63539723/aic-aicc-bic-formula-in-r-for-glm  
  # see also vcdExtra
  loglik <- logLik(model)
  dev <- -2 * as.numeric(loglik)   # LR_chisq
  
  # satd <- loglik + dev / 2
  # LR_chisq <- -2 * (loglik - satd)
  
  df <- model$df.residual
  
  pvalue <- pchisq(dev, df, lower.tail = FALSE)
  
  n   <- attributes(loglik)$nobs
  p   <- attributes(loglik)$df
  
  AIC  <- dev + 2 * p  # model$aic                 -2 * ll + 2 * par
  AICC <- AIC + (2 * p^2 + 2 * p) / (n - p - 1)
  BIC  <- dev +  p * log(n)                   #    -2 * ll + log(ns) * par
  
  outp <- data.frame(loglik, dev, n, df, pvalue, AIC, AICC, BIC)
  colnames(outp)[1:2] <- c('LogLik','Deviance')
  
  return(invisible(outp))
}





# Anova Table (Type III tests)

ANOVA_TABLE <- function(data, model) {
  
  # the lm function uses only Type I sums of squares
  # to get the ANOVA table & partial eta-squared values, need the SS for each
  # term computed after all other terms are in the model
  
  # the procedure below runs the lm function with each predictor
  # having its turn as the final term in the model, in order to get the partial SS
  
  # the Anova function in the car package is probably more efficient, but
  # it does not provide eta-squared
  # the eta_squared function in the effectsize package provides eta-squared but not the anova table
  
  # source for eta_squared = 
  # https://stats.oarc.ucla.edu/stata/faq/how-can-i-compute-effect-size-in-stata-for-regression/
  
  DV <-  names(model.frame(model))[1]
  
  prednoms <- labels(terms(model))
  
  anova_table <- as.data.frame(anova(model))
  
  anova_table <- anova_table[length(prednoms),]
  
  for (lupe in 1:(length(prednoms)-1)) {
    
    new_last <- prednoms[lupe]
    
    rest_noms <- setdiff(prednoms, new_last)
    
    preds_REORD <- c(rest_noms, new_last) 
    
    form_REORD <- as.formula(paste(DV, paste(preds_REORD, collapse=" + "), sep=" ~ "))
    
    model_REORD <- lm(form_REORD, data=data, model=FALSE, x=FALSE, y=FALSE)
    
    # lm: The terms in the formula will be re-ordered so that main effects come first,
    # followed by the interactions, all second-order, all third-order and so on:
    # to avoid this pass a terms object as the formula (see aov and demo(glm.vr) for an example).
    # model_REORD <- lm(mod_terms, data=data, model=FALSE, x=FALSE, y=FALSE)
    # model_REORD <- lm( terms(model), data=data, model=FALSE, x=FALSE, y=FALSE)
    
    anova_table_REORD <- as.data.frame(anova(model_REORD))
    
    anova_table <- rbind(anova_table, anova_table_REORD[length(prednoms),])
  }
  
  Eta_Squared <- anova_table$'Sum Sq'[1:length(prednoms)] / 
    (anova_table$'Sum Sq'[1:length(prednoms)] + anova_table$'Sum Sq'[nrow(anova_table)])
  
  anova_table$Eta_Squared <- Eta_Squared
  
  # reorder the rows to be the same as the orig prednoms
  # use match() to get the row indices that match the prednoms's order
  anova_table <- anova_table[match(prednoms, rownames(anova_table)), ]
  
  return(invisible(anova_table))
}




diagnostics_plots <- function(model, modeldata, plot_diags_nums) {
  
  # # removing quantile plot requests for lm models (they are for glm only)
  # todrop <- intersect(plot_diags_nums, c(10, 11))
  # if (length(todrop) > 0) {
  #   plot_diags_nums <- setdiff(plot_diags_nums, todrop)
  #   cat('\n These plots (in "plot_diags_nums") are not available for lm models:', todrop, '\n')
  # }
  
  Nplots <- length(plot_diags_nums)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (Nplots == 1) par(mfrow = c(1, 1))
  if (Nplots == 2) par(mfrow = c(2, 1),  oma = c(0,4,0,4) )#, mar = c(2,2,1,1))
  if (Nplots == 3) par(mfrow = c(2, 2))
  if (Nplots == 4) par(mfrow = c(2, 2))
  
  for (lupe in 1:length(plot_diags_nums)) {
    
    # portions adapted from
    # Dunn, Peter K. and Gorden K. Smyth. 2018. Generalized Linear Models With Examples in R (Springer).
    
    if (is.element(1, plot_diags_nums)) {
      # fitted (predicted) values vs standardized residuals -- a well fit model will display no pattern.
      plot(modeldata$predicted, modeldata$residuals_standardized,  
           xlab='Predicted (Fitted) Values', ylab='Standardized Residuals',
           main='Standardized Residuals')
      abline(0, 0)
      lowessFit <- lowess(x=modeldata$predicted, y = modeldata$residuals_standardized, 
                          f = 2/3, iter = 3, delta = 0.01 * diff(range(modeldata$predicted)))
      lines(lowessFit,col='red')
    }
    
    if (is.element(2, plot_diags_nums)) {
      # fitted (predicted) values vs studentized residuals
      # https://online.stat.psu.edu/stat462/node/247/#:~:text=Note%20that%20the%20only%20difference,square%20error%20based%20on%20the
      # a studentized residual is the dimensionless ratio resulting from the division 
      # of a residual by an estimate of its standard deviation, both expressed in 
      # the same units. It is a form of a Student's t-statistic, with the estimate 
      # of error varying between points
      plot(modeldata$predicted, modeldata$residuals_studentized,  
           xlab='Predicted (Fitted) Values', ylab='Studentized Residuals',
           main='Studentized Residuals')  #  vs Predicted Values
      abline(0, 0)
      lowessFit <- lowess(x=modeldata$predicted, y = modeldata$residuals_studentized, 
                          f = 2/3, iter = 3, delta = 0.01 * diff(range(modeldata$predicted)))
      lines(lowessFit,col='red')
    }
    
    if (is.element(3, plot_diags_nums)) {
      # histogram of standardized residuals
      hist(modeldata$residuals_standardized, col=4,
           xlab='Standardized Residuals', ylab='Frequency',
           main='Standardized Residuals')  #  Frequencies
    } 
    
    if (is.element(4, plot_diags_nums)) {
      # normal Q-Q plot of standardized residuals
      qqnorm(modeldata$residuals_standardized, pch = 1,
             xlab='Theoretical Values (z)', ylab='Observed Values (z)',
             main='Q-Q Plot: Std. Residuals')
      qqline(modeldata$residuals_standardized, col = 'red', lwd = 2)
    }
    
    # Standard raw residuals arent used in GLM modeling because they dont always make sense.
    # plot(density(resid(model, type='response')))
    
    if (is.element(5, plot_diags_nums)) {
      # density of Pearson residuals
      plot(density(modeldata$residuals_pearson),
           xlab = 'Pearson Residuals',
           main = 'Pearson Residuals')
    }
    
    if (is.element(6, plot_diags_nums)) {
      # density of standardized Pearson residuals -- to ensure the residuals have constant variance
      plot(density(modeldata$residuals_pearson_std),
           xlab = 'Standardized Pearson Residuals',
           main = 'Std. Pearson Residuals')
    }  
    
    if (is.element(7, plot_diags_nums)) {
      # density of deviance residuals, usually preferable to Pearson residuals
      plot(density(modeldata$residuals_deviance),
           xlab = 'Deviance Residuals',
           main = 'Deviance Residuals')
    }
    
    if (is.element(8, plot_diags_nums)) {
      # density of standardized deviance residuals
      plot(density(modeldata$residuals_deviance_std),
           xlab = 'Standardized Deviance Residuals',
           main = 'Std. Deviance Residuals')   
    }
    
    if (is.element(9, plot_diags_nums)) {
      # densities of y and yhat
      # The predicted values from a model will appear similar to y if the model is 
      # well-fit. We hope to see that the distributions are about the same for each model.
      # DV_name <- names(attr(model$terms,"dataClasses")[1])
      DV_name <- names(modeldata)[1]
      
      denplotdat_y     <- density(model$y)
      denplotdat_ypred <- density(modeldata$predicted)
      ymax <- max( c(max(denplotdat_y$y)), max(denplotdat_ypred$y))
      plot(denplotdat_y,
           xlab = DV_name,
           ylim = c(0, ymax),
           main = paste('Original & Predicted ', DV_name, sep=''))
      lines(denplotdat_ypred, col='red')
      if (Nplots < 3)
        legend("topright", legend=c("Orig.", "Pred."), col=c("black", "red"), lty=1, cex=0.8)
    }
    
    if (is.element(10, plot_diags_nums)) {
      # Predicted (Fitted) Values & quantile residuals
      modeldata$residuals_quantile[is.infinite(modeldata$residuals_quantile)] <- NA
      plot(modeldata$predicted, modeldata$residuals_quantile, col='gray',
           xlab = 'Predicted Outcome',
           ylab = 'Quantile Residuals',
           main = 'Quantile Residuals')
      lines(loess.smooth(modeldata$predicted, 
                         modeldata$residuals_quantile), col="red", lty=1, lwd=)
    }
    
    if (is.element(11, plot_diags_nums)) {
      # Q-Q plot with quantile residuals to determine if our chosen distribution makes sense.
      modeldata$residuals_quantile[is.infinite(modeldata$residuals_quantile)] <- NA
      qqnorm(modeldata$residuals_quantile,
             main = 'Q-Q plot: Quantile Residuals')
      qqline(modeldata$residuals_quantile, col = 'red', lwd = 2)
    }
    
    if (is.element(12, plot_diags_nums)) {
      # Cooks distances, to find cases with high leverage   bar plot-like vertical lines
      # same as   plot(model, which = 5, sub.caption='')
      # find the 3 cases with the highest distances
      top3 = order(-modeldata$cooks_distance)[1:3]
      top3_value = modeldata$cooks_distance[top3[3]]
      plot(modeldata$cooks_distance, type='h',
           xlab = 'Case Number',
           ylab = 'Cook\'s Distance',
           main = 'Cook\'s Distances')
      # abline(h=top3_value)
      # labs <- ifelse(modeldata$cooks_distance >= top3_value, row.names(modeldata), "")
      # text(modeldata$cooks_distance, labs, pos=3, col='blue', cex=.8)
    }
    
    # We can further specify the specific points with high leverage using a benchmark 
    # of 2 times the mean of Cooks distance:
    # cooksd_m0 <- cooks.distance(model)  
    # length(cooksd_m0[cooksd_m0 > mean(cooksd_m0) * 2])
    
    # Influence
    # We can also look for influential observations; neither model has any. 
    # which(influence.measures(model)$is.inf[,'cook.d'] ) 
    
    if (is.element(13, plot_diags_nums)) {
      # leverage   bar plot-like vertical lines
      # same as   plot(hatvalues(model), type = 'h')
      # find the 3 cases with the highest leverage values
      # top3 = order(-modeldata$leverage)[1:3]
      # top3_value = modeldata$leverage[top3[3]]
      plot(modeldata$leverage, type = 'h',
           xlab = 'Case Number',
           ylab = 'Leverage',
           main = 'Leverage')
      # labs <- ifelse(modeldata$leverage >= top3_value, row.names(modeldata), "")
      # text(modeldata$leverage, labs, pos=3, col='blue', cex=.8)
    }
    
    if (is.element(14, plot_diags_nums)) {
      # standardized residuals vs leverage
      # same as    plot(model, which = 5, sub.caption='') 
      #calculate leverage for each observation in the model
      plot(modeldata$leverage, modeldata$residuals_standardized,  
           xlab='Leverage', ylab='Standardized Residuals',
           main='Leverage & Std. Residuals')
      lowessFit <- stats::lowess(x=modeldata$leverage, y = modeldata$residuals_standardized, 
                                 f = 2/3, iter = 3, delta = 0.01 * diff(range(modeldata$leverage)))
      lines(lowessFit, col='red')
    }
    
    if (is.element(15, plot_diags_nums)) {
      # Cooks distance vs leverage
      # plot(model, which = 6, sub.caption='') 
      # https://stat.ethz.ch/R-manual/R-patched/library/stats/html/plot.lm.html
      leverage_h = (modeldata$leverage / (1 - modeldata$leverage))
      plot(leverage_h, modeldata$cooks_distance, 
           main = 'Cook\'s Distance vs Leverage*',
           xlab = 'Leverage / (1-Leverage)',
           ylab = 'Cook\'s Distance'
      )
      lowessFit <- stats::lowess(x=leverage_h, y = modeldata$cooks_distance, 
                                 f = 2/3, iter = 3, delta = 0.01 * diff(range(modeldata$leverage)))
      lines(lowessFit, col='red')
    }    
    
    if (is.element(16, plot_diags_nums)) {
      # Scale-Location  sqrt(std residuals) vs fitted -- a well fit model will display no pattern
      plot(modeldata$predicted, sqrt(abs(modeldata$residuals_standardized)),  
           xlab='Predicted (Fitted) Values', ylab='Sqrt. of Standardized Residuals',
           main='Scale-Location')
      abline(0, 0)
      lowessFit <- lowess(x=modeldata$predicted, y = sqrt(abs(modeldata$residuals_standardized)), 
                          f = 2/3, iter = 3, delta = 0.01 * diff(range(modeldata$predicted)))
      lines(lowessFit,col='red')
    }
    
  }    
  #       par(oldpar)
  #       par(mfrow = c(1,1))
  #       dev.off()
  #       graphics.off()
  
}



# adapted from:
# statmod package
# Dunn, Peter K. and Gorden K. Smyth. 2018. Generalized Linear Models With Examples in R (Springer).

quantile_residuals <- function(model) {
  
  y_orig <- model$y
  
  y_fitted <- fitted(model)
  
  N <- length(y_orig)
  
  if (model$family$family == 'binomial' | model$family$family == 'quasibinomial') {
    
    if(!is.null(model$prior.weights)) { n <- model$prior.weights } else { n <- rep(1,N)}
    
    y <- n * y_orig
    mins <- pbinom(y - 1, n, y_fitted)
    maxs <- pbinom(y, n, y_fitted)
  }
  
  if (model$family$family == 'poisson' | model$family$family == 'quasipoisson') {
    
    mins <- ppois(y_orig - 1, y_fitted)
    maxs <- ppois(y_orig,     y_fitted)
  }
  
  if (grepl("Negative Binomial", model$family$family, fixed = TRUE)) {
    
    if (is.null(model$theta)) {size <- model$call$family[[2]]} else {size <- model$theta}
    
    p <- size / (y_fitted + size)
    mins <- ifelse(y_orig > 0, pbeta(p, size, pmax(y_orig, 1)), 0)
    maxs <- pbeta(p, size, y_orig + 1)
  }
  
  u <- runif(n = N, min = mins, max = maxs)
  
  output <- qnorm(u)
  
  return(invisible(output))
}





# portions adapted from 
# https://stackoverflow.com/questions/38109501/how-does-predict-lm-compute-confidence-interval-and-prediction-interval/38110406#38110406

predict_boc <- function(model, modeldata, newdata, CI_level = 95, 
                        bootstrap=FALSE, N_sims=100, kind, family) {
  
  # # computing yhat manually, rather than using built-in predict  NOT USING  prob = factors
  # Xp <- model.matrix(delete.response(terms(model)), newdata)
  # b <- coef(model)
  # yhat <- c(Xp %*% b)  # c() reshape the single-column matrix to a vector
  
  if (kind == 'ZINFL') {
    yhat_count <- predict(model, newdata=newdata, type = "count")
    # zero values are predicted:
    yhat_zero  <- predict(model, newdata=newdata, type = "zero", at = 0)  
    # yhat_zero  <- predict(model, newdata=newdata, type = "prob")  
  } else if (kind == 'HURDLE') {
    yhat_count <- predict(model, newdata=newdata, type = "count")
    # zero values are predicted:
    yhat_zero  <- predict(model, newdata=newdata, type = "prob")[,1]
    # yhat_zero_0  <- predict(model, newdata=newdata, type = "prob", at = 0:4)
    # print(head(yhat_zero))
    # print(head(yhat_zero_0))
  } else {
    yhatR <- predict(model, newdata=newdata, type="response", se.fit = TRUE)
    yhat <- yhatR$fit
    yhat_se.fit <- yhatR$se.fit
  }
  
  
  # CIs
  
  # not doing bootstrapping for MODERATED models because newdata would require values
  # for product terms & lm & glm sometimes alter the var names for interaction terms, too messy
  if (bootstrap & kind == 'MODERATED') {
    message('\nBootstrapped CIs for a MODERATED regression model was requested but it is not')
    message('available in this version of the package. Regular CI values will be provided instead.\n')
    bootstrap <- FALSE
  }  
  
  
  if (!bootstrap) {
    
    if  (kind == 'OLS' | kind == 'MODERATED' | kind == 'LOGISTIC' | 
         kind == 'POISSON' | kind == 'NEGBIN') {
      
      # # computing se.fit manually, rather than using built-in predict  NOT USING  prob = factors
      # var.fit <- rowSums((Xp %*% vcov(model)) * Xp)  # point-wise variance for predicted mean
      # se.fit <- sqrt(var.fit)
      
      zforCI <- qnorm((1 + CI_level * .01) / 2)
      
      ci_lb <- yhat - zforCI * yhat_se.fit
      
      ci_ub <- yhat + zforCI * yhat_se.fit
      
      resmat <- cbind(yhat, ci_lb, ci_ub, yhat_se.fit)
      
      # # not needed when  type="response" on the above predict command
      # if (kind == 'LOGISTIC')
      # # convert to probabilities
      #   resmat <- 1/(1 + exp(-resmat))
      # 
      # if (kind == 'COUNT')
      # # exponential function - inverse of log-link
      #   resmat <- exp(resmat)
      
      # resmat <- cbind(resmat, yhat_se.fit)
    }
    
    if (kind == 'ZINFL' | kind == 'HURDLE') {
      
      ci_lb_count <- ci_ub_count <- ci_lb_zero <- ci_ub_zero <- rep(NA, length(yhat_count))
      
      resmat <- cbind(yhat_count, ci_lb_count, ci_ub_count, yhat_zero, ci_lb_zero, ci_ub_zero)
    }
  }
  
  
  if (bootstrap) {
    
    # simulating parameters using the model -- NOT USING
    # simulate a range of plausible models based on distribution of parameters from model fit 
    # bootparams <- mvrnorm(n=N_sims, mu=coef(model), Sigma=vcov(model))
    
    # real bootstrap parameters
    predicted_values <- predicted_values_count <- predicted_values_zero <- c()
    
    if (kind == 'OLS') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- lm(formula(model), data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }
    
    if (kind == 'LOGISTIC' | kind == 'POISSON') {   
      
      for(i in 1:N_sims) {
        
        bootsamp <- model$data[sample(rownames(model$data), replace=TRUE),]
        
        modelboot <- glm(model$formula, family=family, data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }
    
    if (kind == 'NEGBIN') {   
      
      for(i in 1:N_sims) {
        
        bootsamp <- model$data[sample(rownames(model$data), replace=TRUE),]
        
        modelboot <- MASS::glm.nb(model$formula, 
                                  family=MASS::negative.binomial(1, link="log"), data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }
    
    if (kind == 'ZINFL') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- pscl::zeroinfl(model$formula, dist=family, data=bootsamp)
        
        predicted_values_count <- cbind(predicted_values_count, 
                                        predict(modelboot, newdata=newdata, type='count'))
        
        predicted_values_zero  <- cbind(predicted_values_zero,  
                                        predict(modelboot, newdata=newdata, type='zero'))
      }
    }
    
    if (kind == 'HURDLE') {
      
      # dist	 character specification of count model family.
      # The count model is typically a truncated Poisson or negative binomial regression 
      # (with log link). The geometric distribution is a special case of the 
      # negative binomial with size parameter equal to 1. For modeling the hurdle, 
      # either a binomial model can be employed or a censored count distribution. 
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- pscl::hurdle(model$formula, dist=family, data=bootsamp)
        
        predicted_values_count <- cbind(predicted_values_count, 
                                        predict(modelboot, newdata=newdata, type='count'))
        
        predicted_values_zero  <- cbind(predicted_values_zero,  
                                        predict(modelboot, newdata=newdata, type='zero'))
      }
    }
    
    CI_lo_cut <- (1 - CI_level * .01) / 2
    
    CI_hi_cut <- 1 - (1 - CI_level * .01) / 2
    
    if (kind == 'OLS' | kind == 'LOGISTIC' | kind == 'POISSON' | kind == 'NEGBIN') {
      
      ci_lb <- apply(predicted_values, 1, quantile, probs=CI_lo_cut)
      ci_ub <- apply(predicted_values, 1, quantile, probs=CI_hi_cut)
      
      resmat <- cbind(yhat, ci_lb, ci_ub)
    }
    
    if (kind == 'ZINFL' | kind == 'HURDLE') {
      
      ci_lb_count <- apply(predicted_values_count, 1, quantile, probs=CI_lo_cut)
      ci_ub_count <- apply(predicted_values_count, 1, quantile, probs=CI_hi_cut)
      
      ci_lb_zero  <- apply(predicted_values_zero,  1, quantile, probs=CI_lo_cut)
      ci_ub_zero  <- apply(predicted_values_zero,  1, quantile, probs=CI_hi_cut)
      
      resmat <- cbind(yhat_count, ci_lb_count, ci_ub_count, yhat_zero, ci_lb_zero, ci_ub_zero)
    }
  }
  
  # if (!bootstrap) {
  #    
  #   # comparison with built-in predict
  #   yhatR <- predict(model, newdata=newdata, type="response", se.fit = TRUE)
  #   print(sum(round(c(resmat[,'yhat']   - yhatR$fit),5)))
  #   print(sum(round(c(resmat[,'se.fit'] - yhatR$se.fit),5)))
  # 
  #   # # comparison with ggeffects
  #   # pgg <- ggeffects::predict_response(model, terms = list(AGE = testdata$AGE))
  #   # # print(pgg, n = Inf)
  #   # print(sum(round(c(resmat[,'yhat']   - pgg$predicted),5)))
  #   # print(sum(round(c(resmat[,'ci_lb']  - pgg$conf.low),5)))
  #   # print(sum(round(c(resmat[,'ci_ub']  - pgg$conf.high),5)))
  #   # print(sum(round(c(resmat[,'se.fit'] - pgg$std.error),5)))
  # }
  
  return(invisible(resmat))
}





# rounds numeric columns in a matrix
# numeric columns named 'p' or 'plevel' are rounded to round_p places
# numeric columns not named 'p' are rounded to round_non_p places

round_boc <- function(donnes, round_non_p = 2, round_p = 5) {
  
  # identify the numeric columns
  #	numers <- apply(donnes, 2, is.numeric)  # does not work consistently 
  for (lupec in 1:ncol(donnes)) {
    
    if (is.numeric(donnes[,lupec]) == TRUE) 
      
      if (colnames(donnes)[lupec] == 'p' | colnames(donnes)[lupec] == 'plevel' | 
          colnames(donnes)[lupec] == 'Pr(>|t|)' | colnames(donnes)[lupec] == 'Pr(>F)' )  {
        donnes[,lupec] <- round(donnes[,lupec],round_p)
      } else {
        donnes[,lupec] <- round(donnes[,lupec],round_non_p)				
      }		
    # if (is.numeric(donnes[,lupec]) == FALSE) numers[lupec] = 'FALSE'		
    # if (colnames(donnes)[lupec] == 'p') numers[lupec] = 'FALSE'		
  }
  
  # # set the p column to FALSE
  # numers_not_p <- !names(numers) %in% "p"
  
  #	donnes[,numers_not_p] = round(donnes[,numers_not_p],round_non_p) 
  
  #	if (any(colnames(donnes) == 'p'))  donnes[,'p'] = round(donnes[,'p'],round_p) 
  
  return(invisible(donnes))
}










# Collinearity Diagnostics

Collinearity <- function(modelmat, verbose=FALSE) {
  
  
  # Variance Inflation Factors 
  
  cormat <- cor(modelmat[,-1,drop=FALSE])  # ignoring the intercept column
  
  smc <- 1 - (1 / diag(solve(cormat)))
  
  VIF <- 1 / (1 - smc)
  
  Tolerance <- 1 / VIF
  
  VIFtol <- cbind(Tolerance, VIF)
  
  # Eigenvalues, Condition Indices, Condition Number and Condition Index, & Variance Decomposition Proportions (as in SPSS)
  
  # adapted from https://stat.ethz.ch/pipermail/r-help/2003-July/036756.html
  
  modelmat <- scale(modelmat, center=FALSE) / sqrt(nrow(modelmat) - 1)
  
  svd.modelmat <- svd(modelmat)
  
  singular.values <- svd.modelmat$d
  
  evals <- singular.values **2
  
  condition.indices <- max(svd.modelmat$d) / svd.modelmat$d
  
  phi <- sweep(svd.modelmat$v**2, 2, svd.modelmat$d**2, "/")
  
  VP <- t(sweep(phi, 1, rowSums(phi), "/"))
  
  colnames(VP) <- colnames((modelmat))
  rownames(VP) <- 1:nrow(VP)
  
  
  CondInd <- data.frame(1:length(evals), cbind(evals, condition.indices, VP))	
  colnames(CondInd)[1:4] <- c('Dimension','Eigenvalue','Condition Index', 'Intercept')
  
  if (verbose) {
    message('\n\nCollinearity Diagnostics:\n')
    
    #		message('\n\nTolerance and Variance Inflation Factors:\n')
    print(round(VIFtol, 3), print.gap=4)
    message('\nThe mean Variance Inflation Factor = ', round(mean(VIFtol[,2]),3))
    message('\nMulticollinearity is said to exist when VIF values are > 5, or when Tolerance values')
    message('are < 0.1. Multicollinearity may be biasing a regression model, and there is reason')
    message('for serious concern, when VIF values are greater than 10.\n\n')
    
    print(round(CondInd,3), print.gap=4)
    message('\nThe coefficients for the Intercept and for the predictors are the Variance Proportions.')
    
    message('\nEigenvalues that differ greatly in size suggest multicollinearity. Multicollinearity is')
    message('said to exist when the largest condition index is between 10 and 30, and evidence for')
    message('the problem is strong when the largest condition index is > 30. Multicollinearity is')
    message('also said to exist when the Variance Proportions (in the same row) for two variables are')
    message('above .80 and when the corresponding condition index for a dimension is higher than 10 to 30.\n\n')
  }	
  
  simple.Model <- list(VIFtol=VIFtol, CondInd=CondInd)
  
  return(invisible(simple.Model))
  
  
}

