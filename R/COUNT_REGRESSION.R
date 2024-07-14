

COUNT_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL,
                     family = 'poisson',
                     offset = NULL,
                     plot_type = 'residuals',
                     verbose=TRUE ) {
  
  
  # 2 options in R for Negative Binomial = MASS::glm.nb, & MASS::negative.binomial
  # MASS::negative.binomial replicates SPSS "negative binomial with log link"
  # MASS::glm.nb involves a theta value, and it can be done in SPSS using the Custom 
  # options, but I am not providing this option at this point
  
  # The p-values and S.E.s match those from SPSS when the SPSS
  # parameter estimation method is set to 'Hybrid' and the Scale Parameter Method
  # is set to 'Pearson Chi-square'
  
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if (verbose) {
    message('\nModel Information:')
    message('\nDependent Variable: ', DV)
    message('\nProbability Distribution: ', family)
    if (!is.null(offset)) message('\nThe offset variable is: ', offset)
  }
  
  if (!is.null(forced))  {
    donnes <- data[,c(DV,forced,offset)]	
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (!is.null(hierarchical))  {
    donnes <- data[,c(DV,unlist(hierarchical),offset)]
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
  
  if (is.null(forced) & is.null(hierarchical) ) {
    message('\n\nThe data for the analyses were not properly specified. Expect errors.\n')
    message('\nThe forced & hierarchical arguments were NULL (i.e, not provided).\n')
    message('\nThere is no way of determining what analyses should be conducted.\n')
  }
  
  # gathering all of the predictor variable names
  allIVnoms <- -9999
  if (!is.null(forced))        allIVnoms <- c(allIVnoms, forced)
  if (!is.null(hierarchical))  allIVnoms <- c(allIVnoms, unlist(hierarchical) )
  allIVnoms <- allIVnoms[-1]
  
  donnesRED <- donnes[c(DV,allIVnoms,offset)]  # a version of donnes that contains only the variables in the analyses
  
  
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
  
  
  ### Check contrasts for binary predictors (R creates contrasts automatically
  ### for categorical predictor, with default 0/1 coding)
  # with(couple.df, {
  #   print(contrasts(gender))
  #   print(contrasts(infidelity))
  # })
  # 
  
  if (verbose) {
    message('\nDV frequencies:')
    print(table(donnes[,DV]))
    
    # histogram with a normal curve
    # histogram <- hist(donnes[,DV], breaks=20, col="red", xlab=DV, main="Histogram") 
    
    message('\n\n\nBlock 0: Beginning Block')
  }
  
  # NULL model
  if (is.null(offset))
    formNULL <- as.formula(paste(DV, 1, sep=" ~ "))
  if (!is.null(offset))
    formNULL <- as.formula( paste(paste(DV, 1, sep=" ~ "), 
                                  paste(" + offset(", offset, ")")) )
  
  if (family != 'negbin')
    modelNULL <- glm(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE, 
                     family = family)   # (link="log"))
  
  if (family == 'negbin') {
    # modelNULL <- MASS::glm.nb(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE)
    modelNULL <- glm(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE, 
                     family = MASS::negative.binomial(1, link="log"))
  }
  
  modelNULLsum <- summary(modelNULL)
  null_dev <- data.frame(modelNULL$deviance); colnames(null_dev) <- 'Deviance'
  
  if (verbose) {
    message('\n\nModel Summary:\n')
    print(data.frame(round(null_dev,3)), row.names = FALSE)
  }
  
  # exponentiated coefficients
  exp_B <- exp(modelNULL$coefficients)
  exp_B_CIs <- exp(confint.default(modelNULL))
  
  modelNULLsum$coefs <- cbind(modelNULLsum$coefficients, exp_B, matrix(exp_B_CIs, nrow=1))
  colnames(modelNULLsum$coefs) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
  
  if (verbose) {
    message('\n\nVariables in the Equation -- NULL model:\n')
    modelNULLsum$coefs <- round_boc(modelNULLsum$coefs, round_non_p = 3, round_p = 6) 
    print(round_boc(modelNULLsum$coefs,3), print.gap=4)
  }
  
  prevModel <- modelNULL
  
  if (!is.null(forced)) hierarchical <- list(forced)
  
  for (lupe in 1:length(hierarchical)) {
    
    # keeping info on previous model in lupe for Rsq change stats
    if (lupe > 1) prevModel <- modelMAIN
    
    if (verbose) message('\n\n\nBlock ', lupe)	
    
    if (lupe==1)  preds <- unlist(hierarchical[1])
    
    if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
    
    donnesH <- donnes[,c(DV,preds,offset)]
    
    if (is.null(offset))
      formMAIN <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
    
    if (!is.null(offset))
      formMAIN <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                    paste(" + offset(",  offset, ")")) )
    
    if (family != 'negbin')
      modelMAIN <- glm(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
    
    if (family == 'negbin') {
      modelMAIN <- MASS::glm.nb(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
      # modelMAIN <- glm(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
      #                  family = MASS::negative.binomial(1, link="log")) 
    }
    
    modelMAINsum <- summary(modelMAIN, correlation=TRUE)
    
    # exponentiated coefficients
    exp_B <- exp(modelMAIN$coefficients)
    exp_B_CIs <- exp(confint.default(modelMAIN))
    
    # Goodness of Fit
    X2_Pearson <- sum(residuals(modelMAIN, type = "pearson")^2)
    
    # https://stackoverflow.com/questions/63539723/aic-aicc-bic-formula-in-r-for-glm  
    loglik <- logLik(modelMAIN)
    n   <- attributes(loglik)$nobs
    p   <- attributes(loglik)$df
    dev <- -2*as.numeric(loglik)
    AIC  <- dev + 2 * p  # modelMAIN$aic
    AICC <- AIC + (2 * p^2 + 2 * p)/(n - p - 1)
    BIC  <- dev +  p * log(n)
    
    df <- modelMAIN$df.residual
    
    if (verbose) {
      message('\n\nGoodness of Fit:')
      cat('\n    Deviance: ', round(modelMAIN$deviance,3), '   df =', df, 
          '    value/df =', round((modelMAIN$deviance/df),3)  )
      # message('\n    Scaled Deviance: ')  
      cat('\n\n    Pearson Chi-Square: ', round(X2_Pearson,3), '   df =', df, 
          '    value/df =', round((X2_Pearson/df),3)  )
      # message('\n    Scaled Pearson Chi-Square:')  
      cat('\n\n    Log Likelihood: ', round(loglik,3))  
      cat('\n\n    Akaike Information Criterion (AIC): ', round(AIC,3))  
      cat('\n\n    Finite Sample Corrected AIC (AICC): ', round(AICC,3))  
      cat('\n\n    Bayesian Information Criterion (BIC): ', round(BIC,3))  
      # message('\n    Consistent AIC (CAIC): ')  
    }
    
    # Overdispersion
    
    # Kabacoff p 327: As with logistic regression, overdispersion is suggested 
    # if the ratio of the residual deviance to the residual degrees of freedom 
    # is much larger than 1.
    
    # Overdispersion test based on the model deviance
    dispersion_ratio_dev <- modelMAIN$deviance / modelMAIN$df.residual
    p_deviance <-  1 - pchisq(modelMAIN$deviance, modelMAIN$df.residual)
    if (verbose) {
      message('\n\n\nOverdispersion test based on the model deviance:')
      cat('\n    Dispersion Ratio: ', round(dispersion_ratio_dev,3),
          '    Statistic = ', round(modelMAIN$deviance,3),
          '    p = ', round(p_deviance,5))
    }
    
    # Overdispersion test based on the Pearson Chi-Square
    dispersion_ratio_X2 <- X2_Pearson / modelMAIN$df.residual
    p_X2 <-  1 - pchisq(X2_Pearson, modelMAIN$df.residual)
    if (verbose) {
      message('\n\nOverdispersion test based on the Pearson Chi-Square:')
      cat('\n    Dispersion Ratio: ', round(dispersion_ratio_X2,3),
          '    Statistic = ', round(X2_Pearson,3),
          '    p = ', round(p_X2,5))
    }
    
    # Overdispersion test based on the just the DV
    # adapted from qcc.overdispersion.test
    N <- length(donnes[,DV])
    var_observed <- var(donnes[,DV])
    var_theoretical <- mean(donnes[,DV])
    D <- (var_observed * (N - 1)) / var_theoretical
    ratio_obs.theo <- var_observed / var_theoretical
    pvalue <- 1 - pchisq(D, N - 1)
    if (verbose) {
      message('\n\nOverdispersion test based on just the DV:')
      cat('\n    Obs.Variance / Theor.Variance = ', round(ratio_obs.theo,3), 
          '    Statistic = ', round(D,3), 
          '    p = ', round(pvalue,5))
    }
    
    # Model Effect Sizes
    model_dev <- modelMAIN$deviance 
    null_dev <- modelMAIN$null.deviance 
    model_N <- length(modelMAIN$fitted.values)
    Rsq_HS <-  1 - model_dev / null_dev
    Rsq_CS <- 1- exp (-(null_dev - model_dev) / model_N)
    Rsq_NG <- Rsq_CS / (1 - ( exp(-(null_dev / model_N))))
    model_sumES <- data.frame(Rsq_CS, Rsq_NG, Rsq_HS)
    colnames(model_sumES) <- c('Rsq. Cox & Snell', 'Rsq. Nagelkerke', 'Rsq. Hosmer & Lemeshow')
    if (verbose) {
      message('\n\n\nModel Effect Sizes:\n')
      print(round(model_sumES,3), print.gap=4, row.names = FALSE)
    }
    
    # Omnibus Test
    # compares the fitted model against the intercept-only model
    if (verbose) {
      message('\nOmnibus Test:\n')
      print(anova(prevModel, modelMAIN, test="Chisq"))
    }
    
    # Parameter Estimates
    modelMAINsum$coefs <- cbind(modelMAINsum$coefficients, exp_B, exp_B_CIs)
    colnames(modelMAINsum$coefs) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    
    if (verbose) {
      message('\n\nParameter Estimates:\n')
      modelMAINsum$coefs <- round_boc(modelMAINsum$coefs, round_non_p = 3, round_p = 6) 
      print(round_boc(modelMAINsum$coefs,3), print.gap=4)
    }
    
    # likelihood ratio tests
    if (length(preds) > 1) {
      LR_tests <- c()
      for (lupepreds in 1:length(preds)) {
        
        pred_final <- preds[lupepreds]   
        pred_others <- preds[! preds %in% c(pred_final)]
        
        if (is.null(offset))
          formM1 <- as.formula(paste(DV, paste(pred_others, collapse=" + "), sep=" ~ "))
        if (!is.null(offset))
          formM1 <- as.formula( paste(paste(DV, paste(pred_others, collapse=" + "), sep=" ~ "), 
                                      paste(" + offset(",  offset, ")")) )
        
        if (family != 'negbin')
          M1 <- glm(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (family == 'negbin') {
          # M1 <- MASS::glm.nb(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
          M1 <- glm(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log")) 
        }
        
        if (is.null(offset))
          formM2 <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
        if (!is.null(offset))
          formM2 <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                      paste(" + offset(",  offset, ")")) )
        
        if (family != 'negbin')
          M2 <- glm(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (family == 'negbin') {
          # M2 <- MASS::glm.nb(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
          M2 <- glm(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log"))
        }
        
        aa <- anova(M1, M2, test="Chisq")   # same as from  print(lmtest::lrtest(M1, M2))
        
        LR_tests <- rbind(LR_tests, c(aa$Deviance[2], aa$Df[2], aa$"Pr(>Chi)"[2]))
        #   if (family != 'negbin')
        #     LR_tests <- rbind(LR_tests, c(aa$Deviance[2], aa$Df[2], aa$"Pr(>Chi)"[2]))
        #   if (family == 'negbin')
        #     LR_tests <- rbind(LR_tests, c(aa$"LR stat."[2], aa$"   df"[2], aa$"Pr(Chi)"[2]))
      }
      rownames(LR_tests) <- preds
      colnames(LR_tests) <- c('Likelihood Ratio Chi-Square', 'df', 'p')
      
      
      if (verbose) {
        message('\n\nTests of Model Effects (Type III):\n')
        print(round_boc(LR_tests,3), print.gap=5)
      }
    }
    
    if (verbose) {
      message("\n\nRao's score test (Lagrange Multiplier test):\n")
      # useful to determine whether the Poisson model is appropriate for the data
      print(suppressWarnings(anova(modelMAIN, test = 'Rao')))
      
      message('\n\nCorrelations of Parameter Estimates:\n')
      print(round_boc(modelMAINsum$correlation,3), print.gap=4)
    }
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
  
  
  class(output) <- "COUNT_REGRESSION"
  
  return(invisible(output))
  
}


