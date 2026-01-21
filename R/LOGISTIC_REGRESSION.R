


LOGISTIC_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL, formula=NULL,
                                 ref_category = NULL,
                                 family = 'binomial',
                                 CI_level = 95,
                                 MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                                     thin = 1, burnin = 1000, 
                                                     HDI_plot_est_type = 'standardized'),
                                 plot_type = 'residuals',
                                 verbose=TRUE ) {
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if (verbose) {
    message('\nLogistic Regression:')
    message('\nDependent Variable: ', DV)
    message('\nProbability Distribution: ', family)
  }
  
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
    if (is.null(formula)) 
     formule <- as.formula(paste(DV, paste(forced, collapse=" + "), sep=" ~ "))
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
  allIVnoms <- -9999
  if (!is.null(forced))        allIVnoms <- c(allIVnoms, forced)
  if (!is.null(hierarchical))  allIVnoms <- c(allIVnoms, unlist(hierarchical) )
  allIVnoms <- allIVnoms[-1]
  
  donnes <- donnes[c(DV,allIVnoms)]  # a version of donnes that contains only the variables in the analyses
  
  
  # convert DV to a factor, if it isn't already one
  
  # glm: A typical predictor has the form response ~ terms where response is 
  # the (numeric) response vector and terms is a series of terms which specifies 
  # a linear predictor for response. 
  # For binomial and quasibinomial families the response can also be specified as a 
  # factor (when the first level denotes failure and all others success) or as a 
  # two-column matrix with the columns giving the numbers of successes and failures
  
  # if DV is character, convert it to a factor
    if (is.character(donnes[,DV])) {
      message('\nThe DV is a character variable. It will be converted into a factor.\n')
      donnes[,DV] <- as.factor(donnes[,DV])
    }
  
  # if DV is a factor, confirm that it has 2 levels
  if (is.factor(donnes[,DV])) {
    if (verbose) {
      if (length(levels(donnes[,DV])) == 2) 
        message('\n\nThe DV in data is a factor with two levels.')
      if (length(levels(donnes[,DV])) > 2) 
        message('\n\nThe DV in data is a factor with more than 2 levels. Expect errors.')
    }
  }
  
  # if DV is not a factor, convert it to a factor with 2 levels
  if (!is.factor(donnes[,DV])) {
    donnes[,DV] <- factor(donnes[,DV])
    if (verbose) {
      if (length(levels(donnes[,DV])) == 2) 
        message('\n\nThe DV has been converted to a factor with two levels.')
      if (length(levels(donnes[,DV])) > 2) {
        message('\n\nThe DV has been converted to a factor,')
        message('\nbut there are > 2 levels. Expect errors.')
      }
    }		
  }
  
  if (!is.null(ref_category)) donnes[,DV] <- relevel(donnes[,DV], ref=ref_category)
  
  if (is.null(ref_category))  ref_category <- levels(donnes[,DV])[1]   
  
  if (verbose) {
    message('\nDV frequencies:')
    print(table(donnes[,DV]))
    message('\nThe reference (baseline) category for the DV is: ', ref_category)
  }
  
  
  # clean up IV factors, contrasts
  IV_type_cleanup(donnes, allIVnoms) 
  
  
  # predictor descriptives
  if (verbose)  descriptives(donnes, allIVnoms)
  
  
  
  if (verbose) message('\n\n\nBlock 0: Beginning Block')
  
  # NULL model
  formNULL <- as.formula(paste(DV, 1, sep=" ~ "))
  modelNULL <- glm(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE, family = family)
  modelNULLsum <- summary(modelNULL)
  null_dev <- data.frame(modelNULL$deviance); colnames(null_dev) <- 'Deviance'
  
  if (verbose) {
    message('\n\nModel Summary:\n')
    print(data.frame(round(null_dev,3)), row.names = FALSE)
  }
  
  # Classification Table -- NULL model
  probs_NULL <- fitted(modelNULL)
  classif_tab_NULL <- cbind( c(0,0), table(modelNULL$y, probs_NULL > .5) )
  colnames(classif_tab_NULL) <- rownames(classif_tab_NULL) <- levels(donnes[,DV])
  classif_tab_NULL <- cbind(classif_tab_NULL, c(0,100))
  colnames(classif_tab_NULL)[3] <- 'Percent_Correct'
  names(dimnames(classif_tab_NULL)) <- c('Observed', 'Predicted')
  
  if (verbose) {
    message('\n\nClassification Table for the NULL model:\n')
    print(classif_tab_NULL, print.gap=4)
    message('\nOverall percent correct = ', round((classif_tab_NULL[2,2] / sum(classif_tab_NULL[,2])*100),1))
  }
  
  # odds ratio
  OR <- exp(modelNULL$coefficients)
  OR_CIs <- exp(confint.default(modelNULL))
  
  modelNULLsum$coefs <- cbind(modelNULLsum$coefficients, OR, matrix(OR_CIs, nrow=1))
  colnames(modelNULLsum$coefs) <- c('B', 'SE', 'z', 'p', 'OR', 'OR ci_lb', 'OR ci_ub')
  
  if (verbose) {
    message('\n\nVariables in the Equation -- NULL model:\n')
    modelNULLsum$coefs <- round_boc(modelNULLsum$coefs, round_non_p = 3, round_p = 6) 
    print(round_boc(modelNULLsum$coefs,3), print.gap=4)
  }
  
  prevModel <- modelNULL
  
  if (!is.null(forced)) hierarchical <- list(forced)
  
  for (lupe in 1:length(hierarchical)) {
    
    # keeping info on previous model in lupe for Rsq change stats
    if (lupe > 1)  prevModel <- model

    if (verbose) message('\n\n\nBlock ', lupe)	
    
    if (lupe == 1)  preds <- unlist(hierarchical[1])
    
    if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
    
    # noms_list <- list(DV=DV, IV=preds)
    
    if (is.null(formula))  
      formule <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ ")) 

    model <- glm(formule, data = donnes[,c(DV,preds)], model=TRUE, x=TRUE, y=TRUE, family = family)
    modelsum <- summary(model, correlation=TRUE)
    
    # odds ratios
    OR <- exp(model$coefficients)
    OR_CIs <- exp(confint.default(model))
    
    # Omnibus Tests of Model Coefficients
    # a <- anova(modelNULL, model, test="Chisq")
    if (verbose) {
      message('\n\nOmnibus Tests of Model Coefficients:\n')
      print(anova(prevModel, model, test="Chisq"))
    }
    
    # Overdispersion
    
    # Kabacoff p 327: As with logistic regression, overdispersion is suggested 
    # if the ratio of the residual deviance to the residual degrees of freedom 
    # is much larger than 1.
    
    # # Overdispersion test based on the model deviance
    # dispersion_ratio_dev <- model$deviance / model$df.residual
    # p_deviance <-  1 - pchisq(model$deviance, model$df.residual)
    # message('\n\n\nOverdispersion test based on the model deviance:')
    # cat('\n    Dispersion Ratio: ', round(dispersion_ratio_dev,3),
    #     '    Statistic = ', round(model$deviance,3),
    #     '    p = ', round(p_deviance,5))
    
    # Overdispersion test based on the Pearson Chi-Square
    X2_Pearson <- sum(residuals(model, type = "pearson")^2)
    dispersion_ratio_X2 <- X2_Pearson / model$df.residual
    p_X2 <- 1 - pchisq(X2_Pearson, model$df.residual)
    if (verbose) {
      message('\n\nOverdispersion test:')
      cat('\n    Dispersion Ratio: ', round(dispersion_ratio_X2,3),
          '    Statistic = ', round(X2_Pearson,3),
          '    p = ', round(p_X2,5))
    }
    
    # Model Summary & Effect Sizes
    model_dev <- model$deviance 
    null_dev <- model$null.deviance 
    model_N <- length(model$fitted.values)
    Rsq_HS <-  1 - model_dev / null_dev
    Rsq_CS <- 1- exp (-(null_dev - model_dev) / model_N)
    Rsq_NG <- Rsq_CS / (1 - ( exp(-(null_dev / model_N))))
    # model_sumES <- cbind(model_dev, Rsq_CS, Rsq_NG, Rsq_HS)
    model_sumES <- data.frame(model_dev, Rsq_CS, Rsq_NG, Rsq_HS)
    colnames(model_sumES) <- c('Deviance', 'Rsq. Cox & Snell', 'Rsq. Nagelkerke', 'Rsq. Hosmer & Lemeshow')
    
    if (verbose) {
      message('\n\n\nModel Summary & Effect Sizes:\n')
      print(round(model_sumES,3), print.gap=4, row.names = FALSE)
    }
    
    # Classification Table
    probs <- fitted(model)
    classif_tab <- table(model$y, probs > .5) 
    if (ncol(classif_tab) == 1)  classif_tab <- cbind( c(0,0), table(modelNULL$y, probs_NULL > .5) )
    colnames(classif_tab) <- rownames(classif_tab) <- levels(donnes[,DV])
    names(dimnames(classif_tab)) <- c('Observed', 'Predicted')
    
    if (verbose) {
      message('\n\nClassification Table:\n')
      print(classif_tab, print.gap=4)
      message('\nThe cut value for the classifications is .5\n')
    }
    
    modelsum$coefs <- cbind(modelsum$coefficients, OR, OR_CIs)
    colnames(modelsum$coefs) <- c('B', 'SE', 'z', 'p', 'OR', 'OR ci_lb', 'OR ci_ub')
    
    if (verbose) {
      message('\nVariables in the Equation:\n')
      modelsum$coefs <- round_boc(modelsum$coefs, round_non_p = 3, round_p = 6) 
      print(round_boc(modelsum$coefs,3), print.gap=4)
    }
    
    MCMC_outp <- list(Bayes_HDIs = NULL, chains = NULL, autocorrels = NULL, SDs = NULL)
    
    if (MCMC_options$MCMC) {
      
      if (family == 'quasibinomial') {
        message("\n\nFamily = 'quasibinomial' analyses are currently not possible for")
        message("the MCMC analyses. family = 'binomial' will therefore be used instead.\n")
      }

      MCMC_outp <- Bayes_HDIs_reg_preds(formule, data = donnes[,c(DV,preds)], CI_level,  
                                        MCMC_options, model_type = 'LOGISTIC')  # , noms_list,
        
      if (verbose) {
        message('\n\nBayesian MCMC chains:')
        message('\n    Number iterations (samples): ', MCMC_options$Nsamples)
        message('\n    Thinning interval:  ', MCMC_options$thin)
        message('\n    Burn-in period:  ', MCMC_options$burnin)
        message('\n    Final chain length (# of samples): ', nrow(MCMC_outp$chains))
        message('\n    The priors were the default rstanarm::stan_glm function priors')
        
        message('\n\nAutocorrelations for the MCMC chains:\n')
        print(round(MCMC_outp$autocorrels,3), print.gap=4)
        
        message('\n\nBayesian HDIs for the raw coefficients & the corresponding odds ratios:\n')
        print(round_boc(MCMC_outp$Bayes_HDIs, round_non_p = 2, round_p = 5), print.gap=4)
      }
    }
    
    # likelihood ratio tests
    if (length(preds) > 1) {
      LR_tests <- c()
      for (lupepreds in 1:length(preds)) {
        
        pred_final <- preds[lupepreds]   
        pred_others <- preds[! preds %in% c(pred_final)]
        formM1 <- as.formula(paste(DV, paste(pred_others, collapse=" + "), sep=" ~ "))
        M1 <- glm(formM1, data = donnes[,c(DV,preds)], model=TRUE, x=TRUE, y=TRUE, family = family)
        formM2 <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
        M2 <- glm(formM2, data = donnes[,c(DV,preds)], model=TRUE, x=TRUE, y=TRUE, family = family)
        # print(lmtest::lrtest(M1, M2))
        aa <- anova(M1, M2, test="Chisq")
        LR_tests <- rbind(LR_tests, c(aa$Deviance[2], aa$Df[2], aa$"Pr(>Chi)"[2]))
      }
      rownames(LR_tests) <- preds 
      colnames(LR_tests) <- c('Likelihood Ratio Chi-Square', 'df', 'p')
      
      if (verbose) {
        message('\n\nTests of Model Effects (Type III):\n')
        print(round_boc(LR_tests,3), print.gap=5)
      }
    }
    
    if (verbose) {
      message('\n\nCorrelations of Parameter Estimates:\n')
      print(round_boc(modelsum$correlation,3), print.gap=4)
    }
    
    modeldata <- data.frame(model$y, model$x[,2:ncol(model$x)])
    colnames(modeldata) <- c(DV,colnames(model$x)[2:ncol(model$x)])
    
    modeldata$predicted <- model$fitted.values
    
    # casewise diagnostics
    modeldata$residuals_raw <-          resid(model, type='response')
    modeldata$residuals_pearson <-      resid(model, type='pearson')
    modeldata$residuals_pearson_std <-  rstandard(model, type='pearson')
    modeldata$residuals_deviance <-     resid(model, type='deviance')
    modeldata$residuals_deviance_std <- rstandard(model, type='deviance')
    modeldata$residuals_quantile <-     quantile_residuals(model)
    
    modeldata$cooks_distance <- cooks.distance(model)
    modeldata$dfbeta <- dfbeta(model)
    modeldata$dffit <- dffits(model)
    modeldata$leverage <- hatvalues(model)
    modeldata$covariance_ratios <- covratio(model)
    
    collin_diags <- Collinearity(model.matrix(model), verbose=FALSE)
  }
  
  
  if (verbose) {			
    message('\n\nPredictor variable correlation matrix:\n')
    print(round(cor(modeldata[,c(DV,names(model$coefficients)[-1])]),3), print.gap=4)		
    
    message('\n\nCollinearity Diagnostics:')
    VIFtol_messages(collin_diags$VIFtol[,2])
    CondInd_messages(collin_diags$CondInd)
  }     
  
  
  # add, if any, factor variables in data to modeldata
  # (because lm changes names & types and the original variables are needed for PLOTMODEL)
  # factor_variables <- names(model$model[sapply(model$model, is.factor)])
  factor_variables <- list_xlevels <- names(model$xlevels)
  if (!is.null(factor_variables))  modeldata[factor_variables] <- donnes[,factor_variables]
  
  
  output <- list(model=model, modelsum=modelsum,
                 modeldata=modeldata, coefs = modelsum$coefs,
                 collin_diags = collin_diags, family=family,
                 chain_dat = MCMC_outp$chains,
                 Bayes_HDIs = MCMC_outp$Bayes_HDIs)  #noms_list = noms_list)
                
  
  if (plot_type == 'residuals') 
    
    diagnostics_plots(model=model, modeldata=modeldata, plot_diags_nums=c(8, 9, 10, 11))
  
  if (plot_type == 'diagnostics')
    
    diagnostics_plots(model=model, modeldata=modeldata, plot_diags_nums=c(12, 13, 14, 15))
  
  if (plot_type == 'Bayes_HDI' & MCMC_options$MCMC & !is.null(MCMC_outp$chain)) 
    
    Bayes_HDI_plot(chain_dat = MCMC_outp$chains, 
                   Bayes_HDIs = MCMC_outp$Bayes_HDIs,
                   CI_level = CI_level, # SDs = xx$SDs[c(preds)],
                   HDI_plot_est_type = MCMC_options$HDI_plot_est_type)
  
  
  class(output) <- "LOGISTIC_REGRESSION"
  
  return(invisible(output))
}


