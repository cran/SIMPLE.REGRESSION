

# forced=NULL; hierarchical=NULL;
# plot_type = 'residuals';
# CI_level = 95;
# MCMC = FALSE;
# Nsamples = 10000;
# verbose=TRUE
# 
# data=data_Green_Salkind_2014; DV='injury';
#                forced = c('quads','gluts','abdoms','arms','grip')



OLS_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL,
                            plot_type = 'residuals', 
                            CI_level = 95,
                            MCMC = FALSE,
                            Nsamples = 10000,
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
  donnesMOD <- donnes[c(DV,allIVnoms)]  
  
  # converting character variables to factors
  for (lupe in 1:length(allIVnoms)) {
    if (is.character(donnesMOD[,allIVnoms[lupe]])) {
      message('\n', allIVnoms[lupe], ' is a character variable. It will be converted into a factor.\n')
      donnesMOD[,allIVnoms[lupe]] <- as.factor(donnesMOD[,allIVnoms[lupe]])
    }
  }

  # check if all IVs are either numeric or factors
  IVclasses <- sapply(donnesMOD[,allIVnoms], class)
  if (!all(IVclasses == 'numeric' | IVclasses == 'integer' | IVclasses == 'factor')) {
    message('\n The following variables are not numeric, integers, or factors. Expect errors.')
    names(
    IVclasses[IVclasses != 'numeric' &  IVclasses != 'integer' & IVclasses != 'factor' ])
  }

  # check if DV is numeric
  if (!is.numeric(donnesMOD[,DV])) message('\n DV is not numeric. Expect errors.')
  
  
  # descriptives
  if (verbose) {
    # descriptives for numeric variables
    donnesNUM <- donnesMOD[sapply(donnesMOD,is.numeric)] # selecting only numeric variables
    if (ncol(donnesNUM) != 0) {
      minmax <- t(apply(donnesNUM, 2, range))
      descs <- data.frame(Mean=colMeans(donnesNUM), SD=apply(donnesNUM, 2, sd), 
                          Min=minmax[,1], Max=minmax[,2]) 
      message('\n\nDescriptive statistics for the numeric variables:\n')
      print(round(descs,2), print.gap=4)
    }
    
    # frequencies for factors
    donnesFAC <- donnesMOD[sapply(donnesMOD,is.factor)] # selecting only factor variables
    if (ncol(donnesFAC) != 0) {
      message('\n\nCategory frequencies for the factor variables:\n')
      print(apply((donnesFAC), 2, table))
    }
    
    rm(donnesNUM, donnesFAC)
  }
  
  
  # test if factor variables have contrasts; if not, add contrasts with names
  IV_factors <- IVclasses[IVclasses == 'factor' ]
  if (length(IV_factors) > 0) {
    for (lupe in 1:length(IV_factors)) {
      
      if (is.null(attr( donnesMOD[,names(IV_factors[lupe])], "contrasts" ))) {
        
        fac_levels <- levels(donnesMOD[,names(IV_factors[lupe])])
        
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
          better_contrast_noms(contrasts(donnesMOD[,names(IV_factors[lupe])]))
        
        # apply the contrast matrix to the factor variable
        contrasts(donnesMOD[,names(IV_factors[lupe])]) <- custom_contrasts
      }
    }
  }
    
  # display the contrasts for factors
  if (verbose & length(IV_factors) > 0) {
    for (lupe in 1:length(IV_factors)) {
      message('\nThe contrasts for ', names(IV_factors[lupe]), ' are:\n')
      print(contrasts(donnesMOD[,names(IV_factors[lupe])]), print.gap=4)
    }
  }
  
  

  # non hierarchical regression (without interactions or blocks)
  if (!is.null(forced))  {
    
    formMAIN <- as.formula(paste(DV, paste(forced, collapse=" + "), sep=" ~ "))
    
    modelMAIN <- lm(formMAIN, data=donnesMOD, model=TRUE, x=TRUE, y=TRUE)
    
    modelMAINsum <- summary(modelMAIN)
    
    RsqMAIN <- modelMAINsum$r.squared
    
    # Using the classical R^2 test statistic for (linear) regression designs, 
    # this function computes the corresponding Bayes factor test
    # the Bayes factor (against the intercept-only null)
    Nvars <- ncol(donnesMOD) - 1
    BF_mod_Rsq <- linearReg.R2stat(R2 = RsqMAIN, N = nrow(donnesMOD), p = Nvars, 
                                   simple = TRUE)
    
    # Anova Table (Type III tests)
    anova_table <- ANOVA_TABLE(data=donnesMOD, model=modelMAIN) 
    
    # creating a new version of the raw data, from the lm output, so that can 
    # provide stats for dummy codes
    # modeldata is needed for moderation analyses (next)
    modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
    colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
    
    mainRcoefs <- 
      PARTIAL_COR(data=cor(modeldata), Y=DV, X=names(modelMAIN$coefficients)[-1], 
                  Ncases=nrow(donnesMOD), verbose=FALSE)
    
    mainRcoefs <- cbind(mainRcoefs$betas, mainRcoefs$Rx_y, 
                        mainRcoefs$R_partials, mainRcoefs$R_semipartials)
    colnames(mainRcoefs) <- c('beta','r','partial.r','semipartial.r')
    
    
    modeldata$predicted <- modelMAIN$fitted.values
    
    # casewise diagnostics
    modeldata$residuals <- resid(modelMAIN)
    modeldata$residuals_standardized <- rstandard(modelMAIN)
    modeldata$residuals_studentized <- rstudent(modelMAIN)
    modeldata$cooks_distance <- cooks.distance(modelMAIN)
    modeldata$dfbeta <- dfbeta(modelMAIN)
    modeldata$dffit <- dffits(modelMAIN)
    modeldata$leverage <- hatvalues(modelMAIN)
    modeldata$covariance_ratios <- covratio(modelMAIN)
    
    collin_diags <- Collinearity(model.matrix(modelMAIN), verbose=FALSE)
    
    if (MCMC) {
      
      mod_lmBF <- lmBF(formMAIN, data = donnesMOD) 
      
      # Bayes factors - model   Against denominator: Intercept only
      BF_mod_h1h0 <- unlist(extractBF(mod_lmBF)[1])
      BF_mod_h0h1 <- 1 / BF_mod_h1h0
      
      chains <- posterior(mod_lmBF, iterations = Nsamples, progress = FALSE)
      
      Nests <- ncol(chains) - 2
      
      # parameters estimates & CIs
      ests <- colMeans(chains[,1:Nests])
      quant_size <- (100 - CI_level) / 2
      quant_lb <- quant_size * .01
      quant_ub <- (100 - quant_size) * .01
      ests_ci_lb <- apply(chains[,1:Nests],2,quantile,probs=quant_lb)
      ests_ci_ub <- apply(chains[,1:Nests],2,quantile,probs=quant_ub)
      
      # https://cran.r-project.org/web/packages/BayesFactor/vignettes/manual.html#regression
      # The results are quite similar (to OLS coefs), apart from the intercept. This is due to the 
      # Bayesian model centering the covariates before analysis, so the mu parameter 
      # is the mean of  rather than the expected value of the response variable 
      # when all uncentered covariates are equal to 0.
      
      # standardized estimates & CIs
      SD_preds <- apply(donnesMOD,2,sd)
      ests_z <- c(NA, (ests[-1] * SD_preds[-1]) / SD_preds[1])
      ests_z_ci_lb <- c(NA, (ests_ci_lb[-1] * SD_preds[-1]) / SD_preds[1])
      ests_z_ci_ub <- c(NA, (ests_ci_ub[-1] * SD_preds[-1]) / SD_preds[1])
      
      # Bayes factors - for the estimates, from just t & N
      t_values <- modelMAINsum$coefficients[,3]
      BF_ests_h1h0 <- unlist(sapply(t_values, ttest.tstat, n1 = nrow(donnesMOD), 
                                    simple=TRUE))
      BF_ests_h0h1 <- 1 / BF_ests_h1h0
      
      Bayes_ests <- cbind(ests, ests_ci_lb, ests_ci_ub, 
                          ests_z, ests_z_ci_lb, ests_z_ci_ub,
                          BF_ests_h1h0, BF_ests_h0h1)
      rownames(Bayes_ests)[1] <- 'Intercept'
    }      
    
    
    if (verbose) {	
      
      message('\n\nThe DV is: ', DV)
      message('\nThe IVs are: ', paste(forced, collapse=', ') )
      
      message('\n\nmultiple R = ', round(sqrt(RsqMAIN),3),  
              '   multiple R-squared = ', round(RsqMAIN,3),
              '   adjusted R-squared = ', round(modelMAINsum$adj.r.squared,3))
      
      Fstats <- modelMAINsum$fstatistic
      pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
      
      message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2],  
              '   df_denom = ', Fstats[3], '   p-value = ', round(pvalue,6), '\n')
      
      BF_interps(BF_10 = BF_mod_Rsq, BF_01 = 1 / BF_mod_Rsq) 
      
      message('\n\nAnova Table (Type III tests):\n')
      print(round_boc(anova_table,3), print.gap=4)
      
      modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, confint(modelMAIN))
      message('\n\nModel Coefficients:\n')
      print(round_boc(modelMAINsum$coefficients,3), print.gap=4)
      
      message('\n\nBeta, r, partial correlations, & semi-partial correlations:\n')
      print(round(mainRcoefs,3), print.gap=4)	
      
      if (MCMC) {
        
        # message('\n\nBayes Factor for model vs. null: ', round(BF_mod_h1h0,2))
        # message('\nBayes Factor for null vs. model: ', round(BF_mod_h0h1,2))
        # message('\nBayes Factor for model based on Rsq.: ', round(BF_mod_Rsq,2))
        
        message('\n\nBayesian Raw and Standardized Coefficients and Bayes Factors:\n')
        
        colnames(Bayes_ests) <- cbind('b', 'b_ci_lb', 'b_ci_ub', 
                                      'Beta', 'Beta_ci_lb', 'Beta_ci_ub',
                                      'BF_10', 'BF_01')
        
        print(round_boc(Bayes_ests, round_non_p = 2, round_p = 5), print.gap=4)
      }
      
      message('\n\nCorrelation matrix:\n')
      # if (!is.null(forced)) print(round(cor(modeldata[,c(DV,forced)]),3), print.gap=4)		
      print(round(cor(modeldata[,c(DV,names(modelMAIN$coefficients)[-1])]),3), print.gap=4)		
      
      message('\n\nCollinearity Diagnostics:\n')
      print(round(collin_diags$VIFtol, 3), print.gap=4)
      message('\nThe mean Variance Inflation Factor = ', round(mean(collin_diags$VIFtol[,2]),3))
      message('\nMulticollinearity is said to exist when VIF values are > 5, or when Tolerance values')
      message('are < 0.1. Multicollinearity may be biasing a regression model, and there is reason')
      message('for serious concern, when VIF values are greater than 10.\n\n')
      print(round(collin_diags$CondInd,3), print.gap=4)
      message('\nThe coefficients for the Intercept and for the predictors are the Variance Proportions.')
      message('\nEigenvalues that differ greatly in size suggest multicollinearity. Multicollinearity is')
      message('said to exist when the largest condition index is between 10 and 30, and evidence for')
      message('the problem is strong when the largest condition index is > 30. Multicollinearity is')
      message('also said to exist when the Variance Proportions (in the same row) for two variables are')
      message('above .80 and when the corresponding condition index for a dimension is higher than 10 to 30.\n\n')
    }
    
    
    # add, if any, factor variables in data to modeldata 
    # (because lm changes names & types and the original variables are needed for PLOTMODEL)
    factor_variables <- names(modelMAIN$model[sapply(modelMAIN$model, is.factor)])
    if (!is.null(factor_variables))  modeldata[factor_variables] <- donnesMOD[,factor_variables]
    
    output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum, 
                   anova_table=anova_table, mainRcoefs=mainRcoefs, 
                   modeldata=modeldata, collin_diags=collin_diags, family='OLS')
    
    class(output) <- "OLS_REGRESSION"
    
  }
  
  
  
  # hierarchical regression (blocks/steps)
  if (!is.null(hierarchical)) {
    
    for (lupe in 1:length(hierarchical)) {
      
      # keeping info on previous model in lupe for Rsq change stats
      if (lupe > 1) {
        prevModel <- modelMAIN
        prevRsq   <- RsqMAIN
        if (MCMC) prevmod_lmBF <- mod_lmBF
      }
      
      message('\n\nStep ', lupe)	
      
      if (lupe==1)  preds <- unlist(hierarchical[1])
      
      if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
      
      donnesH <- donnesMOD[,c(DV,preds)]
      
      formMAIN <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
      
      modelMAIN <- lm(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE)
      
      modelMAINsum <- summary(modelMAIN)
      
      RsqMAIN <- modelMAINsum$r.squared
      
      # Using the classical R^2 test statistic for (linear) regression designs, 
      # this function computes the corresponding Bayes factor test
      # the Bayes factor (against the intercept-only null)
      Nvars <- ncol(donnesH) - 1
      BF_mod_Rsq <- linearReg.R2stat(R2 = RsqMAIN, N = nrow(donnesH), p = Nvars, 
                                     simple = TRUE)
      
      # Anova Table (Type III tests)
      anova_table <- ANOVA_TABLE(data=donnesH, model=modelMAIN) 
      
      # creating a new version of the raw data, from the lm output, 
      # so that can provide stats for dummy codes
      modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
      colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
      
      mainRcoefs <- 
        PARTIAL_COR(data=cor(modeldata), Y=DV, X=names(modelMAIN$coefficients)[-1], 
                    Ncases=nrow(donnesH), verbose=FALSE)
      
      mainRcoefs <- cbind(mainRcoefs$betas, mainRcoefs$Rx_y, 
                          mainRcoefs$R_partials, mainRcoefs$R_semipartials)
      colnames(mainRcoefs) <- c('beta','r','partial.r','semipartial.r')
      
      
      if (MCMC) {
        
        mod_lmBF <- lmBF(formMAIN, data = donnesH)
        
        chains <- posterior(mod_lmBF, iterations = Nsamples, progress = FALSE)
        
        Nests <- ncol(chains) - 2
        
        # parameters estimates & CIs
        ests <- colMeans(chains[,1:Nests])
        quant_size <- (100 - CI_level) / 2
        quant_lb <- quant_size * .01
        quant_ub <- (100 - quant_size) * .01
        ests_ci_lb <- apply(chains[,1:Nests],2,quantile,probs=quant_lb)
        ests_ci_ub <- apply(chains[,1:Nests],2,quantile,probs=quant_ub)
        
        # https://cran.r-project.org/web/packages/BayesFactor/vignettes/manual.html#regression
        # The results are quite similar, apart from the intercept. This is due to the 
        # Bayesian model centering the covariates before analysis, so the mu parameter 
        # is the mean of  rather than the expected value of the response variable 
        # when all uncentered covariates are equal to 0.
        
        # standardized estimates & CIs
        SD_preds <- apply(donnesH,2,sd)
        ests_z <- c(NA, (ests[-1] * SD_preds[-1]) / SD_preds[1])
        ests_z_ci_lb <- c(NA, (ests_ci_lb[-1] * SD_preds[-1]) / SD_preds[1])
        ests_z_ci_ub <- c(NA, (ests_ci_ub[-1] * SD_preds[-1]) / SD_preds[1])
        
        # Bayes factors - for the estimates, from just t & N
        t_values <- modelMAINsum$coefficients[,3]
        BF_ests_h1h0 <- unlist(sapply(t_values, ttest.tstat, n1 = nrow(donnesH), simple=TRUE))
        BF_ests_h0h1 <- 1 / BF_ests_h1h0
        
        Bayes_ests <- cbind(ests, ests_ci_lb, ests_ci_ub, 
                            ests_z, ests_z_ci_lb, ests_z_ci_ub,
                            BF_ests_h1h0, BF_ests_h0h1)
        rownames(Bayes_ests)[1] <- 'Intercept'
      }
      
      modeldata$predicted <- modelMAIN$fitted.values
      
      # casewise diagnostics
      modeldata$residuals <- resid(modelMAIN)
      modeldata$residuals_standardized <- rstandard(modelMAIN)
      modeldata$residuals_studentized <- rstudent(modelMAIN)
      modeldata$cooks_distance <- cooks.distance(modelMAIN)
      modeldata$dfbeta <- dfbeta(modelMAIN)
      modeldata$dffit <- dffits(modelMAIN)
      modeldata$leverage <- hatvalues(modelMAIN)
      modeldata$covariance_ratios <- covratio(modelMAIN)
      
      collin_diags <- Collinearity(model.matrix(modelMAIN), verbose=FALSE)
      
      if (verbose) {	
        
        message('\nThe DV is: ', DV); message('\nThe IVs are: ', paste(preds, collapse=', ') )
        
        message('\n\nmultiple R = ', round(sqrt(RsqMAIN),3), 
                '   multiple R-squared = ', round(RsqMAIN,3),
                '   adjusted R-squared = ', round(modelMAINsum$adj.r.squared,3))
        
        Fstats <- modelMAINsum$fstatistic
        pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
        
        message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2],  
                '   df_denom = ', Fstats[3], '   p-value = ', round(pvalue,6), '\n')
        
        BF_interps(BF_10 = BF_mod_Rsq, BF_01 = (1 / BF_mod_Rsq) )
        
        # if (MCMC) {
        #   
        #   message('\n\nBayes Factor for model vs. null: ', round(BF_mod_h1h0,2))
        #   message('\nBayes Factor for null vs. model: ', round(BF_mod_h0h1,2))
        #   message('\nBayes Factor for model based on Rsq.: ', round(BF_mod_Rsq,2))
        # }
        
        if (lupe > 1) {
          # current vs previous model comparisons	
          Rsqch <- RsqMAIN - prevRsq
          message('\n\nCurrent vs previous model comparison:')
          message('\n    Rsquared change = ', round(Rsqch,3))
          
          fish <- anova(modelMAIN, prevModel)
          message('\n    F = ', round(fish$F[2],2), '   df_num = ', 1,  
                  '   df_denom = ', fish$Res.Df[1], '   p-value = ', 
                  round(fish$'Pr(>F)'[2],6), '\n')
          
          if (MCMC) {
            
            prevmod_lmBF@data <- mod_lmBF@data
            
            BF_change <- mod_lmBF / prevmod_lmBF
            
            BF_change <- as.numeric(unlist(extractBF(BF_change))[1])
            
            BF_interps(BF_M2 = BF_change)
            
          }
        }
        
        message('\n\nAnova Table (Type III tests):\n')
        print(round_boc(anova_table,3), print.gap=4)
        
        message('\nModel Coefficients:\n')
        modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, confint(modelMAIN))
        print(round_boc(modelMAINsum$coefficients,3), print.gap=4)
        
        message('\nBeta, r, partial correlations, & semi-partial correlations:\n')
        print(round(mainRcoefs,3), print.gap=4)
        
        if (MCMC) {
          
          message('\nBayesian Raw and Standardized Coefficients and Bayes Factors:\n')
          
          colnames(Bayes_ests) <- cbind('b', 'b_ci_lb', 'b_ci_ub', 
                                        'Beta', 'Beta_ci_lb', 'Beta_ci_ub',
                                        'BF_10', 'BF_01')
          
          print(round_boc(Bayes_ests, round_non_p = 2, round_p = 5), print.gap=4)
        }
        
      }
    }	
    
    if (verbose) {			
      message('\n\nCorrelation matrix:\n')
      
      print(round(cor(modeldata[,c(DV,names(modelMAIN$coefficients)[-1])]),3), print.gap=4)		
      
      message('\n\nCollinearity Diagnostics:\n')
      print(round(collin_diags$VIFtol, 3), print.gap=4)
      message('\nThe mean Variance Inflation Factor = ', round(mean(collin_diags$VIFtol[,2]),3))
      message('\nMulticollinearity is said to exist when VIF values are > 5, or when Tolerance values')
      message('are < 0.1. Multicollinearity may be biasing a regression model, and there is reason')
      message('for serious concern, when VIF values are greater than 10.\n\n')
      print(round(collin_diags$CondInd,3), print.gap=4)
      message('\nThe coefficients for the Intercept and for the predictors are the Variance Proportions.')
      message('\nEigenvalues that differ greatly in size suggest multicollinearity. Multicollinearity is')
      message('said to exist when the largest condition index is between 10 and 30, and evidence for')
      message('the problem is strong when the largest condition index is > 30. Multicollinearity is')
      message('also said to exist when the Variance Proportions (in the same row) for two variables are')
      message('above .80 and when the corresponding condition index for a dimension is higher than 10 to 30.\n\n')
    }
    
    
    # add, if any, factor variables in data to modeldata 
    # (because lm changes names & types and the original variables are needed for PLOTMODEL)
    factor_variables <- names(modelMAIN$model[sapply(modelMAIN$model, is.factor)])
    if (!is.null(factor_variables))  modeldata[factor_variables] <- donnesMOD[,factor_variables]
    
    
    output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum, 
                   anova_table=anova_table, mainRcoefs=mainRcoefs, modeldata=modeldata, 
                   collin_diags=collin_diags, family='OLS')
    
    class(output) <- "OLS_REGRESSION"
    
  }
  
  
  if (plot_type == 'residuals') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(16,2,3,4))
  
  if (plot_type == 'diagnostics') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(9,12,13,14))
  
  
  return(invisible(output))
  
}


SIMPLE.REGRESSION <- OLS_REGRESSION





