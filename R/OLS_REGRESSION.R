

OLS_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL, formula=NULL, 
                            CI_level = 95,
                            MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                                thin = 1, burnin = 1000, 
                                                HDI_plot_est_type = 'standardized'),
                            plot_type = 'residuals', 
                            verbose=TRUE, ... ) {
  
  args <- as.list(sys.call())
  if ("MOD" %in% names(args)) {
    stop('\nThe package functions have been changed:
          \nplease use MODERATED_REGRESSION for moderated regressions\n\n')
  }
  
  if (as.character(match.call()[[1]]) == "SIMPLE.REGRESSION") {
    warning('\nThe package functions have been changed:
             \nplease use OLS_REGRESSION instead of SIMPLE.REGRESSION for simultaneous and hierarchical regressions
             \nplease use MODERATED_REGRESSION for moderated regressions\n\n', call. = FALSE) 
  }
  
  # "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if (!is.null(formula)) {
    formula_vars <- formula_check(formula=formula, data = data)
    DV <- formula_vars$DV
    forced <- formula_vars$forced
    formule <- formula_vars$formule
    hierarchical <- NULL
  }
  
  if (!is.null(forced))  {
    donnes <- data[,c(DV,forced)]	
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (!is.null(hierarchical))  {
    donnes <- data[,c(DV,unlist(hierarchical))]
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (is.null(forced) & is.null(hierarchical) ) {
    message('\n\nThe data for the analyses were not properly specified.\n')
    message('\nThe forced & hierarchical arguments were both NULL (i.e, not provided). Expect errors.\n')
  }
  
  if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
  
  
  # gathering all of the predictor variable names
  allIVnoms <- c()
  if (!is.null(forced))        allIVnoms <- c(allIVnoms, forced)
  if (!is.null(hierarchical))  allIVnoms <- c(allIVnoms, unlist(hierarchical))
  
  # a version of donnes that contains only the variables in the analyses
  donnes <- donnes[c(DV,allIVnoms)]  
  
  # clean up IV factors, contrasts
  IV_type_cleanup(donnes, allIVnoms) 
    
  # predictor descriptives
  if (verbose)  descriptives(donnes, allIVnoms)

  

  if (!is.null(forced)) hierarchical <- list(forced)
  
  for (lupe in 1:length(hierarchical)) {
    
    # keeping info on previous model in lupe for Rsq change stats
    if (lupe > 1) {
      prevModel <- xx$model
      prevRsq   <- xx$Rsq
      if (MCMC_options$MCMC) prevmod_lmBF <- xx$mod_lmBF
    }
    
    if (lupe == 1)  prevRsq <- prevModel <- prevmod_lmBF <- NULL
    
    message('\n\nBlock ', lupe, ':')	
    
    if (lupe==1)  preds <- unlist(hierarchical[1])
    
    if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
    
    # noms_list <- list(DV=DV, IV=preds)
    
    if (is.null(formula)) 
      formule <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ ")) 
    
    xx <- reg_stats_boc(formule = formule, data = donnes[,c(DV,preds)], 
                        prevRsq = prevRsq, prevModel = prevModel, prevmod_lmBF = prevmod_lmBF,
                        CI_level = CI_level, MCMC_options = MCMC_options) #, noms_list = noms_list)
    
    if (verbose) {	
      
      message('\nThe DV is: ', DV); message('\nThe IVs are: ', paste(preds, collapse=', ') )
      
      resultats(xx)
    }
  }
  
  
  # add, if any, factor variables in data to modeldata 
  # (because lm changes names & types and the original variables are needed for PLOTMODEL)
  factor_variables <- names(xx$model$model[sapply(xx$model$model, is.factor)])
  if (!is.null(factor_variables))  xx$modeldata[factor_variables] <- donnes[,factor_variables]
  
  
  output <- list(model=xx$model, modelsum=xx$modelsum, 
                 anova_table=xx$anova_table, partialRcoefs=xx$partialRcoefs, modeldata=xx$modeldata, 
                 collin_diags=xx$collin_diags, family='OLS',
                 chain_dat = xx$chains, 
                 Bayes_HDIs = xx$Bayes_HDIs) #noms_list = noms_list)
  
  class(output) <- "OLS_REGRESSION"
  
  
  if (plot_type == 'residuals')
    diagnostics_plots(model=xx$model, modeldata=xx$modeldata, plot_diags_nums=c(16,2,3,4))
  
  if (plot_type == 'diagnostics')
    diagnostics_plots(model=xx$model, modeldata=xx$modeldata, plot_diags_nums=c(9,12,13,14))
  
  if (plot_type == 'Bayes_HDI' & MCMC_options$MCMC & !is.null(xx$chains)) 
    Bayes_HDI_plot(chain_dat = xx$chains, 
                   Bayes_HDIs = xx$Bayes_HDIs,
                   CI_level = CI_level, # SDs = xx$SDs[c(preds)],
                   HDI_plot_est_type = MCMC_options$HDI_plot_est_type)
  
  return(invisible(output))
}


SIMPLE.REGRESSION <- OLS_REGRESSION



