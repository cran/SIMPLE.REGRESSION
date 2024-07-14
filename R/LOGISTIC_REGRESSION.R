

LOGISTIC_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL,
                      ref_category = NULL,
                      family = 'binomial',
                      plot_type = 'residuals',
                      verbose=TRUE ) {
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if (verbose) {
    message('\nLogistic Regression:')
    message('\nDependent Variable: ', DV)
    message('\nProbability Distribution: ', family)
  }
  
  if (!is.null(forced))  {
    donnes <- data[,c(DV,forced)]	
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
    formMAIN <- as.formula(paste(DV, paste(forced, collapse=" + "), sep=" ~ "))
  }
  
  if (!is.null(hierarchical))  {
    donnes <- data[,c(DV,unlist(hierarchical))]
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
  
  if (is.null(forced) & is.null(hierarchical)) {
    message('\n\nThe data for the analyses were not properly specified. Expect errors.\n')
    message('\nThe forced & hierarchical arguments were NULL (i.e, not provided).\n')
    message('\nThere is no way of determining what analyses should be conducted.\n')
  }
  
  # gathering all of the predictor variable names
  allIVnoms <- -9999
  if (!is.null(forced))        allIVnoms <- c(allIVnoms, forced)
  if (!is.null(hierarchical))  allIVnoms <- c(allIVnoms, unlist(hierarchical) )
  allIVnoms <- allIVnoms[-1]
  
  donnesRED <- donnes[c(DV,allIVnoms)]  # a version of donnes that contains only the variables in the analyses
  
  
  # descriptives
  if (verbose) {
    # descriptives for numeric variables
    donnesNUM <- donnes[sapply(donnesRED,is.numeric)] # selecting only numeric variables
    if (ncol(donnesNUM) != 0) {
      minmax <- t(apply(donnesNUM, 2, range))
      descs <- data.frame(Mean=colMeans(donnesNUM), SD=apply(donnesNUM, 2, sd), Min=minmax[,1], Max=minmax[,2]) 
      message('\n\nDescriptive statistics for the numeric variables:\n')
      print(round(descs,2), print.gap=4)
    }
    
    # frequencies for factors
    donnesFAC <- donnes[sapply(donnesRED,is.factor)] # selecting only factor variables
    if (ncol(donnesFAC) != 0) {
      message('\n\nCategory frequencies for the factor variables:\n')
      print(apply((donnesFAC), 2, table))
    }
    rm(donnesRED, donnesNUM, donnesFAC)
  }
  
  # convert DV to a factor, if it isn't already one
  
  # glm: A typical predictor has the form response ~ terms where response is 
  #the (numeric) response vector and terms is a series of terms which specifies 
  # a linear predictor for response. 
  # For binomial and quasibinomial families the response can also be specified as a 
  # factor (when the first level denotes failure and all others success) or as a 
  # two-column matrix with the columns giving the numbers of successes and failures
  
  # convert DV to factor if it is non numeric & has 2 levels
  if (is.factor(donnes[,DV])) {
    if (verbose) {
      if (length(levels(donnes[,DV])) == 2) 
        message('\n\nThe DV in data is a factor with two levels.')
      if (length(levels(donnes[,DV])) > 2) 
        message('\n\nThe DV in data has more than 2 levels. The analyses cannot be performed.')
    }
  }
  if (!is.factor(donnes[,DV])) {
    donnes[,DV] <- factor(donnes[,DV])
    if (verbose) {
      if (length(levels(donnes[,DV])) == 2) 
        message('\n\nThe DV in data is non numeric and has been converted to a factor with two levels.')
      if (length(levels(donnes[,DV])) > 2) {
        message('\n\nThe DV in data has been converted to a factor,')
        message('\nbut there are > 2 levels. The analyses cannot be performed.')
      }
    }		
  }
  
  if (verbose) {
    message('\nDV frequencies:')
    print(table(donnes[,DV]))
  }
  
  # specify the DV baseline category
  # ref_category = 'did not graduate'
  # ref_category = 'NO'
  
  if (!is.null(ref_category))  donnes[,DV] <- relevel(donnes[,DV], ref = ref_category)
  
  if ( is.null(ref_category))  ref_category <- levels(donnes[,DV])[1]   
  
  if (verbose) 
    message('\nThe reference (baseline) category for the DV is: ', ref_category)
  
  
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
    if (lupe > 1) {
      prevModel <- modelMAIN
      # prevRsq   <- RsqMAIN
    }
    
    if (verbose) message('\n\n\nBlock ', lupe)	
    
    if (lupe==1)  preds <- unlist(hierarchical[1])
    
    if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
    
    donnesH <- donnes[,c(DV,preds)]
    
    formMAIN <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
    
    modelMAIN <- glm(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
    modelMAINsum <- summary(modelMAIN, correlation=TRUE)
    
    # odds ratios
    OR <- exp(modelMAIN$coefficients)
    OR_CIs <- exp(confint.default(modelMAIN))
    
    # Omnibus Tests of Model Coefficients
    # a <- anova(modelNULL, modelMAIN, test="Chisq")
    if (verbose) {
      message('\n\nOmnibus Tests of Model Coefficients:\n')
      print(anova(prevModel, modelMAIN, test="Chisq"))
    }
    
    # Overdispersion
    
    # Kabacoff p 327: As with logistic regression, overdispersion is suggested 
    # if the ratio of the residual deviance to the residual degrees of freedom 
    # is much larger than 1.
    
    # # Overdispersion test based on the model deviance
    # dispersion_ratio_dev <- modelMAIN$deviance / modelMAIN$df.residual
    # p_deviance <-  1 - pchisq(modelMAIN$deviance, modelMAIN$df.residual)
    # message('\n\n\nOverdispersion test based on the model deviance:')
    # cat('\n    Dispersion Ratio: ', round(dispersion_ratio_dev,3),
    #     '    Statistic = ', round(modelMAIN$deviance,3),
    #     '    p = ', round(p_deviance,5))
    
    # Overdispersion test based on the Pearson Chi-Square
    X2_Pearson <- sum(residuals(modelMAIN, type = "pearson")^2)
    dispersion_ratio_X2 <- X2_Pearson / modelMAIN$df.residual
    p_X2 <- 1 - pchisq(X2_Pearson, modelMAIN$df.residual)
    if (verbose) {
      message('\n\nOverdispersion test:')
      cat('\n    Dispersion Ratio: ', round(dispersion_ratio_X2,3),
          '    Statistic = ', round(X2_Pearson,3),
          '    p = ', round(p_X2,5))
    }
    
    # Model Summary & Effect Sizes
    model_dev <- modelMAIN$deviance 
    null_dev <- modelMAIN$null.deviance 
    model_N <- length(modelMAIN$fitted.values)
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
    probs <- fitted(modelMAIN)
    classif_tab <- table(modelMAIN$y, probs > .5) 
    if (ncol(classif_tab) == 1)  classif_tab <- cbind( c(0,0), table(modelNULL$y, probs_NULL > .5) )
    colnames(classif_tab) <- rownames(classif_tab) <- levels(donnes[,DV])
    names(dimnames(classif_tab)) <- c('Observed', 'Predicted')
    
    if (verbose) {
      message('\n\nClassification Table:\n')
      print(classif_tab, print.gap=4)
      message('\nThe cut value for the classifications is .5\n')
    }
    
    modelMAINsum$coefs <- cbind(modelMAINsum$coefficients, OR, OR_CIs)
    colnames(modelMAINsum$coefs) <- c('B', 'SE', 'z', 'p', 'OR', 'OR ci_lb', 'OR ci_ub')
    
    if (verbose) {
      message('\nVariables in the Equation:\n')
      modelMAINsum$coefs <- round_boc(modelMAINsum$coefs, round_non_p = 3, round_p = 6) 
      print(round_boc(modelMAINsum$coefs,3), print.gap=4)
    }
    
    # likelihood ratio tests
    if (length(preds) > 1) {
      LR_tests <- c()
      for (lupepreds in 1:length(preds)) {
        
        pred_final <- preds[lupepreds]   
        pred_others <- preds[! preds %in% c(pred_final)]
        formM1 <- as.formula(paste(DV, paste(pred_others, collapse=" + "), sep=" ~ "))
        M1 <- glm(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
        formM2 <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
        M2 <- glm(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
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
      print(round_boc(modelMAINsum$correlation,3), print.gap=4)
    }
    
    modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
    colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
    
    modeldata$predicted <- modelMAIN$fitted.values
    
    # casewise diagnostics
    modeldata$residuals_raw <-          resid(modelMAIN, type='response')
    modeldata$residuals_pearson <-      resid(modelMAIN, type='pearson')
    modeldata$residuals_pearson_std <-  rstandard(modelMAIN, type='pearson')
    modeldata$residuals_deviance <-     resid(modelMAIN, type='deviance')
    modeldata$residuals_deviance_std <- rstandard(modelMAIN, type='deviance')
    modeldata$residuals_quantile <-     quantile_residuals(modelMAIN)
    
    modeldata$cooks_distance <- cooks.distance(modelMAIN)
    modeldata$dfbeta <- dfbeta(modelMAIN)
    modeldata$dffit <- dffits(modelMAIN)
    modeldata$leverage <- hatvalues(modelMAIN)
    modeldata$covariance_ratios <- covratio(modelMAIN)
  
    collin_diags <- Collinearity(model.matrix(modelMAIN), verbose=FALSE)
  }
  
  
  if (verbose) {			
    message('\n\nPredictor variable correlation matrix:\n')
    #		print(round(cor(modeldata[,c(DV,preds)]),3), print.gap=4)	
    
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
  
  output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum,
                 modeldata=modeldata, coefs = modelMAINsum$coefs,
                 cormat = modelMAINsum$correlation,
                 collin_diags = collin_diags, family=family)
  
  
  if (plot_type == 'residuals') 
    
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(8, 9, 10, 11))
    
  if (plot_type == 'diagnostics')
    
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(12, 13, 14, 15))
  
 
  class(output) <- "LOGISTIC_REGRESSION"
  
  return(invisible(output))
}


