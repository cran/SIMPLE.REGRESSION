

COUNT_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL,
                              model_type = 'poisson',
                              offset = NULL,
                              plot_type = 'residuals',
                              CI_level = 95,
                              GoF_model_types = TRUE,
                              verbose=TRUE ) {

  # # MCMC = FALSE, Nsamples = 4000,  Sept 2025: the rstanarm package is being removed from CRAN
  
  
  # 2 options in R for Negative Binomial = MASS::glm.nb, & MASS::negative.binomial
  # MASS::negative.binomial replicates SPSS "negative binomial with log link"
  # MASS::glm.nb involves a theta value, and it can be done in SPSS using the Custom 
  # options, but I am not providing this option at this point
  
  # The p-values and S.E.s match those from SPSS when the SPSS
  # parameter estimation method is set to 'Hybrid' and the Scale Parameter Method
  # is set to 'Pearson Chi-square'
  
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  # check of whether the specified model_type is correct/possible
  if (model_type != 'poisson' & model_type != 'quasipoisson' & model_type != 'negbin' & 
      model_type != 'zinfl_poisson' & model_type != 'zinfl_negbin' & 
      model_type != 'hurdle_poisson' & model_type != 'hurdle_negbin') 
    message('\nError: The specified model_type is not one of the possible options.\n')
  
  # create family & kind from model_type
  
  if (grepl('poisson', model_type, fixed=TRUE))       family <- 'poisson'
  if (grepl('quasipoisson', model_type, fixed=TRUE))  family <- 'quasipoisson'
  if (grepl('negbin', model_type, fixed=TRUE))        family <- 'negbin'
  
  if (grepl('zinfl',  model_type, fixed=TRUE))                 kind <- 'ZINFL'
  if (grepl('hurdle', model_type, fixed=TRUE))                 kind <- 'HURDLE'
  if (model_type == 'poisson' | model_type == 'quasipoisson')  kind <- 'POISSON'
  if (model_type == 'negbin')                                  kind <- 'NEGBIN'
  
  if (verbose) {
    message('\nModel Information:')
    message('\nDependent Variable: ', DV)
    message('\nModel type: ', model_type)
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
  
  if (NAflag) cat('\nCases with missing values were found and removed from the data matrix.\n')
  
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
  }
  
  
  if (verbose)  message('\n\n\nBlock 0: Beginning Block')
  
  # NULL model
  if (is.null(offset))
    formNULL <- as.formula(paste(DV, 1, sep=" ~ "))
  
  if (!is.null(offset))
    formNULL <- as.formula( paste(paste(DV, 1, sep=" ~ "), 
                                  paste(" + offset(", offset, ")")) )
  
  if (kind == 'POISSON')
    modelNULL <- glm(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE, 
                     family = family)   # (link="log"))
  
  if (kind == 'NEGBIN') {
    # modelNULL <- MASS::glm.nb(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE)
    modelNULL <- glm(formNULL, data = donnes, model=TRUE, x=TRUE, y=TRUE, 
                     family = MASS::negative.binomial(1, link="log"))
  }
  
  if (kind == 'ZINFL')
    modelNULL <- pscl::zeroinfl(formNULL, data = donnes, 
                                model=TRUE, x=TRUE, y=TRUE, dist = family)
  
  if (kind == 'HURDLE')
    modelNULL <- pscl::hurdle(formNULL, data = donnes, 
                              model=TRUE, x=TRUE, y=TRUE, dist = family)
  
  if (kind == 'ZINFL' | kind == 'HURDLE') 
    modelNULL$deviance <- -2 * modelNULL$loglik
  
  modelNULLsum <- summary(modelNULL)
  
  null_dev <- data.frame(modelNULL$deviance); colnames(null_dev) <- 'Deviance'
  
  if (verbose) {
    message('\n\nModel Summary:\n')
    print(data.frame(round(null_dev,3)), row.names = FALSE)
  }
  
  # exponentiated coefficients
  
  if (kind == 'POISSON' | kind == 'NEGBIN') {
    
    exp_B <- exp(modelNULL$coefficients)
    exp_B_CIs <- exp(confint.default(modelNULL))
    
    modelNULLsum$coefs <- cbind(modelNULLsum$coefficients, exp_B, matrix(exp_B_CIs, nrow=1))
    colnames(modelNULLsum$coefs) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    
    if (verbose) {
      message('\n\nVariables in the Equation -- NULL model:\n')
      modelNULLsum$coefs <- round_boc(modelNULLsum$coefs, round_non_p = 3, round_p = 6) 
      print(round_boc(modelNULLsum$coefs,3), print.gap=4)
    }
  }
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {         # FIXXX?   
    
    exp_B <- exp(unlist(modelNULL$coefficients))
    
    exp_B_CIs <- exp(confint.default(modelNULL))
    
    modelNULLsum$coefs$count <- 
      cbind( matrix(modelNULLsum$coefficients$count[1,],nrow=1), exp_B[1], matrix(exp_B_CIs[1,], nrow=1))
    colnames(modelNULLsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    
    modelNULLsum$coefs$zero <- 
      cbind( matrix(modelNULLsum$coefficients$zero[1,],nrow=1), exp_B[2], matrix(exp_B_CIs[2,], nrow=1))
    colnames(modelNULLsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    
    if (verbose) {
      message('\n\nVariables in the Equation -- NULL model, count portion:\n')
      modelNULLsum$coefs$count <- round_boc(modelNULLsum$coefs$count, round_non_p = 3, round_p = 6) 
      print(round_boc(modelNULLsum$coefs$count,3), print.gap=4)
      
      message('\n\nVariables in the Equation -- NULL model, zero portion:\n')
      modelNULLsum$coefs$zero <- round_boc(modelNULLsum$coefs$zero, round_non_p = 3, round_p = 6) 
      print(round_boc(modelNULLsum$coefs$zero,3), print.gap=4)
    }
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
    
    if (kind == 'POISSON')
      modelMAIN <- glm(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
    
    if (kind == 'NEGBIN') {
      modelMAIN <- MASS::glm.nb(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
      # modelMAIN <- glm(formMAIN, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
      #                  family = MASS::negative.binomial(1, link="log")) 
    }
    
    if (kind == 'ZINFL')
      modelMAIN <- pscl::zeroinfl(formMAIN, data = donnesH, 
                                  model=TRUE, x=TRUE, y=TRUE, dist = family)
    
    if (kind == 'HURDLE')
      modelMAIN <- pscl::hurdle(formMAIN, data = donnes, 
                                model=TRUE, x=TRUE, y=TRUE, dist = family)
    
    if (kind == 'ZINFL' | kind == 'HURDLE')
      modelMAIN$deviance <- -2 * modelMAIN$loglik
    
    modelMAINsum <- summary(modelMAIN, correlation=TRUE)
    
    # exponentiated coefficients, & adding them to modelMAINsum$coefs
    
    if (kind == 'POISSON' | kind == 'NEGBIN') {
      exp_B <- exp(modelMAIN$coefficients)
      exp_B_CIs <- exp(confint.default(modelMAIN))
      
      modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, exp_B, exp_B_CIs)
      colnames(modelMAINsum$coefficients) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    }
    
    if (kind == 'ZINFL' | kind == 'HURDLE') {
      
      exp_B_count <- exp(modelMAIN$coefficients$count)
      exp_B_zero  <- exp(modelMAIN$coefficients$zero)
      
      exp_B_CIs <- exp(confint.default(modelMAIN))
      
      # dropping the row from the count coefs
      modCcoefs <- modelMAINsum$coefficients$count[row.names(modelMAINsum$coefficients$count) != "Log(theta)",]
      modelMAINsum$coefs$count <- 
        cbind(modCcoefs, exp_B_count, exp_B_CIs[grepl("count", rownames(exp_B_CIs)),])
      colnames(modelMAINsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
      
      modelMAINsum$coefs$zero <- 
        cbind(modelMAINsum$coefficients$zero, exp_B_zero, exp_B_CIs[grepl("zero", rownames(exp_B_CIs)),])
      colnames(modelMAINsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    }
    
    
    # Goodness of Fit
    GoFs <- GoF_stats(model = modelMAIN)
    
    if (verbose) {
      message('\n\nGoodness of Fit:')
      
      cat('\n    Log Likelihood: ', round(GoFs$LogLik,3))  
      
      cat('\n\n    Deviance: ', round(GoFs$Deviance,3), '   df =', GoFs$df, 
          '    p value =', round(GoFs$pvalue,5)  )
      
      # cat('\n    Deviance: ', round(GoFs$Deviance,3), '   df =', GoFs$df, 
      #     '    p value =', round((GoFs$Deviance / GoFs$df),3)  )
      # message('\n    Scaled Deviance: ')  
      
      # cat('\n\n    Pearson Chi-Square: ', round(GoFs$X2_Pearson,3), '   df =', GoFs$df, 
      #     '    value/df =', round((GoFs$X2_Pearson / GoFs$df),3)  )
      # message('\n    Scaled Pearson Chi-Square:')  
      
      cat('\n\n    Akaike Information Criterion (AIC): ', round(GoFs$AIC,3))  
      cat('\n\n    Finite Sample Corrected AIC (AICC): ', round(GoFs$AICC,3))  
      cat('\n\n    Bayesian Information Criterion (BIC): ', round(GoFs$BIC,3))  
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
    X2_Pearson <- sum(residuals(modelMAIN, type = "pearson")^2)
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
    # null_dev <- modelMAIN$null.deviance 
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
    
    # Comparisons of the current model against the previous model
    if (verbose) {
      
      if (kind == 'POISSON' | kind == 'NEGBIN') {
        message('\nComparisons of the current model against the previous model:\n')
        print(anova(prevModel, modelMAIN, test="Chisq"))
      }
      
      if (kind == 'ZINFL' |  kind == 'HURDLE') {
        
        # # LRT to compare these models
        # lmtest::lrtest(modelMAIN, prevModel)
        # lmtest::waldtest(modelMAIN, prevModel)
        
        # Compare the current model to a null model without predictors using 
        # chi-squared test on the difference of log likelihoods
        diff_logliks <- as.numeric(2 * (logLik(modelMAIN) - logLik(prevModel)))
        p <- pchisq(diff_logliks, df = length(preds), lower.tail = FALSE) 
        message('\nComparisons of the current model against the previous model:\n')
        message('\n   Difference in the log likelihoods = ', round(diff_logliks,3), 
                '    df = ', length(preds), '     p = ' , round(p,8), '\n\n')
        
        # Which model should we choose?   https://rpubs.com/lucymark2013/817199
        print(pscl::vuong(modelMAIN, prevModel))
      }
    }
    
    # Parameter Estimates
    if (verbose) {
      
      if (kind == 'POISSON' | kind == 'NEGBIN') {
        message('\n\nParameter Estimates:\n')
        modelMAINsum$coefs <- round_boc(modelMAINsum$coefficients, round_non_p = 3, round_p = 6) 
        print(round_boc(modelMAINsum$coefs,3), print.gap=4)
      }
      
      if (kind == 'ZINFL' | kind == 'HURDLE') {
        message('\n\nVariables in the Equation -- count portion of the model:\n')
        modelMAINsum$coefs$count <- round_boc(modelMAINsum$coefs$count, round_non_p = 3, round_p = 6) 
        print(round_boc(modelMAINsum$coefs$count,3), print.gap=4)
        
        message('\n\nVariables in the Equation -- zero portion of the model:\n')
        modelMAINsum$coefs$zero <- round_boc(modelMAINsum$coefs$zero, round_non_p = 3, round_p = 6) 
        print(round_boc(modelMAINsum$coefs$zero,3), print.gap=4)
        
        # Zeileis - Regression Models for Count Data in R, p 19:
        # For the hurdle model, the zero hurdle component describes the probability 
        # of observing a positive count whereas, for the ZINB model, the 
        # zero-inflation component predicts the probability of observing a zero 
        # count from the point mass component.
        
        if (kind == 'ZINFL') {
          message('\n\nThe zero portion of the model is predicting the probability of observing a zero count.')
          message('A positive coefficient (B) for a predictor thus means that as values on a predictor increase,')
          message('the probability of observing a zero value for ', DV, ' increases.')
        }
        if (kind == 'HURDLE') {
          message('\n\nThe zero portion of the model is predicting the probability of observing a non-zero count.') 
          message('A positive coefficient (B) for a predictor thus means that as values on a predictor increase,')
          message('the probability of crossing the hurdle (obtaining a value higher than zero) for ', DV, ' increases.')
        }
      }
    }
    
    # # Sept 2025: the rstanarm package is being removed from CRAN
    # if (MCMC & (kind == 'POISSON' | kind == 'NEGBIN')) {
    #   
    #   if (family == 'poisson') 
    #     MCMC_mod <- rstanarm::stan_glm(formMAIN, data = donnesH, family = family, 
    #                                    refresh = 0, algorithm="sampling", iter = Nsamples)
    #   
    #   if (family == 'quasipoisson') {
    #     message("\n\nFamily = 'quasipoisson' analyses are currently not possible for")
    #     message("the MCMC analyses. family = 'poisson' will therefore be used instead.\n")
    #     MCMC_mod <- rstanarm::stan_glm(formMAIN, data = donnesH, family = "poisson", 
    #                                    refresh = 0, algorithm="sampling", iter = Nsamples)
    #   }
    #   
    #   if (kind == 'NEGBIN') 
    #     MCMC_mod <- rstanarm::stan_glm(formMAIN, data = donnesH, family = "neg_binomial_2", 
    #                                    refresh = 0, algorithm="sampling", iter = Nsamples)
    #   
    #   MCMC_mod_coefs <- cbind(coef(MCMC_mod), 
    #                           rstanarm::posterior_interval(MCMC_mod, prob= CI_level * .01)[1:length(coef(MCMC_mod)),]   )
    #   
    #   MCMC_mod_coefs <- cbind(MCMC_mod_coefs, exp(MCMC_mod_coefs))
    #   
    #   colnames(MCMC_mod_coefs) <- c('B', 'B_ci_lb', 'B_ci_ub',
    #                                 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    #   
    #   if (verbose) {
    #     message('\n\nCoefficients from Bayesian MCMC analyses:\n')
    #     print(round_boc(MCMC_mod_coefs,3), print.gap=4)
    #   }
    # }
    
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
        
        if (kind == 'POISSON')
          M1 <- glm(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (kind == 'NEGBIN') {
          # M1 <- MASS::glm.nb(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
          M1 <- glm(formM1, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log")) 
        }
        
        if (kind == 'ZINFL')
          M1 <- pscl::zeroinfl(formM1, data = donnesH, 
                               model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (kind == 'HURDLE')
          M1 <- pscl::hurdle(formM1, data = donnesH, 
                             model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (is.null(offset))
          formM2 <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
        
        if (!is.null(offset))
          formM2 <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                      paste(" + offset(",  offset, ")")) )
        
        if (kind == 'POISSON')
          M2 <- glm(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (kind == 'NEGBIN') {
          # M2 <- MASS::glm.nb(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE)
          M2 <- glm(formM2, data = donnesH, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log"))
        }
        
        if (kind == 'ZINFL')
          M2 <- pscl::zeroinfl(formM2, data = donnesH, 
                               model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (kind == 'HURDLE')
          M2 <- pscl::hurdle(formM2, data = donnesH, 
                             model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (kind == 'POISSON' | kind == 'NEGBIN') {
          
          aa <- anova(M1, M2, test="Chisq")   # same as from  print(lmtest::lrtest(M1, M2))
          
          LR_tests <- rbind(LR_tests, c(aa$Deviance[2], aa$Df[2], aa$"Pr(>Chi)"[2]))
          #   if (family != 'negbin')
          #     LR_tests <- rbind(LR_tests, c(aa$Deviance[2], aa$Df[2], aa$"Pr(>Chi)"[2]))
          #   if (kind == 'NEGBIN')
          #     LR_tests <- rbind(LR_tests, c(aa$"LR stat."[2], aa$"   df"[2], aa$"Pr(Chi)"[2]))
        }
        
        if (kind == 'ZINFL' | kind == 'HURDLE') {
          
          # Compare the current model to a null model without predictors using 
          # chi-squared test on the difference of log likelihoods
          diff_logliks <- as.numeric(2 * (logLik(M2) - logLik(M1)))
          
          M1_sum <- summary(M1)
          M2_sum <- summary(M2)
          df_comps <- abs(M2_sum$df.residual - M1_sum$df.residual)
          
          p <- pchisq(diff_logliks, df = df_comps, lower.tail = FALSE) 
          # message('\nComparisons of the current model against the previous model:\n')
          # message('\n   Difference in the log likelihoods = ', round(diff_logliks,3), 
          #         '    df = ', df_comps, '     p = ' , round(p,8), '\n\n')
          
          LR_tests <- rbind(LR_tests, c(diff_logliks, df_comps, p))
        }
      }
      rownames(LR_tests) <- preds
      colnames(LR_tests) <- c('Likelihood Ratio Chi-Square', 'df', 'p')
      
      if (verbose) {
        message('\n\nTests of Model Effects (Type III):\n')
        print(round_boc(LR_tests,3), print.gap=5)
      }
    }
    
    if (verbose & (kind == 'POISSON' | kind == 'NEGBIN')) {
      message("\n\nRao's score test (Lagrange Multiplier test):\n")
      # useful to determine whether the Poisson model is appropriate for the data
      print(suppressWarnings(anova(modelMAIN, test = 'Rao')))
      
      message('\n\nCorrelations of Parameter Estimates:\n')
      print(round_boc(modelMAINsum$correlation,3), print.gap=4)
    }
  }
  
  
  # modeldata
  if (kind == 'POISSON')  modeldata <- modelMAIN$data
  
  if (kind == 'NEGBIN') {
    
    if (is.null(offset)) {
      
      modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
      
      colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
    }
    
    if (!is.null(offset)) {
      
      modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)], modelMAIN$offset)
      
      colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)], offset)
    }
  }
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    # modeldata <- data.frame(modelMAIN$y, modelMAIN$x$zero[,2:ncol(modelMAIN$x$zero)])
    # colnames(modeldata) <- c(DV,colnames(modelMAIN$x$zero)[2:ncol(modelMAIN$x$zero)])
    modeldata <- modelMAIN$model
    
    if (!is.null(offset))  colnames(modeldata)[ncol(modeldata)] <- offset
  }
  
  modeldata$predicted <- modelMAIN$fitted.values
  
  
  # modelcoefs <- modelMAINsum$coefs
  modelcoefs <- modelMAINsum$coefficients
  
  
  # variable correlations - using model data with dummy coded factor variables (all numeric)
  if (kind == 'POISSON' | kind == 'NEGBIN') {
    modeldataN <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
    colnames(modeldataN) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
  }
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    modeldataN <- data.frame(modelMAIN$y, modelMAIN$x$zero[,2:ncol(modelMAIN$x$zero)])
    colnames(modeldataN) <- c(DV,colnames(modelMAIN$x$zero)[2:ncol(modelMAIN$x$zero)])
  }
  modeldatacorrels <- cor(modeldataN)
  
  
  # casewise diagnostics
  modeldata$residuals_raw <-          resid(modelMAIN, type='response')
  modeldata$residuals_pearson <-      resid(modelMAIN, type='pearson')
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    modeldata$residuals_pearson_std <-  NULL
    modeldata$residuals_deviance <-     NULL
    modeldata$residuals_deviance_std <- NULL
    modeldata$residuals_quantile <-     NULL
    modeldata$cooks_distance <-         NULL
    modeldata$dfbeta <-                 NULL
    modeldata$dffit <-                  NULL
    modeldata$leverage <-               NULL
    modeldata$covariance_ratios <-      NULL
  } else { 
    modeldata$residuals_pearson_std <-  rstandard(modelMAIN, type='pearson') 
    modeldata$residuals_deviance <-     resid(modelMAIN, type='deviance')
    modeldata$residuals_deviance_std <- rstandard(modelMAIN, type='deviance')
    modeldata$residuals_quantile <-     quantile_residuals(modelMAIN)
    modeldata$cooks_distance <-         cooks.distance(modelMAIN)
    modeldata$dfbeta <-                 dfbeta(modelMAIN)
    modeldata$dffit <-                  dffits(modelMAIN)
    modeldata$leverage <-               hatvalues(modelMAIN)
    modeldata$covariance_ratios <-      covratio(modelMAIN)
  }
  
  collin_diags <- Collinearity(model.matrix(modelMAIN), verbose=FALSE)
  
  if (verbose) {			
    message('\n\nVariable correlation matrix:\n')
    print(round(modeldatacorrels,3), print.gap=4)
    
    #	print(round(cor(modeldata[,c(DV,preds)]),3), print.gap=4)	
    # print(round(cor(modeldata[,c(DV,names(modelMAIN$coefficients)[-1])]),3), print.gap=4)		
    
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
  
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(5, 9))
    
  } else {
    
    if (plot_type == 'residuals') 
      
      diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(8, 9, 10, 11))
    
    if (plot_type == 'diagnostics')
      
      diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(12, 13, 14, 15))
  }
  
  
  if (GoF_model_types) {
    # Goodness of Fit for all possible model types
    
    preds <- unlist(hierarchical[1])
    
    donnesH <- donnes[,c(DV,preds,offset)]
    
    if (is.null(offset))
      formMAIN <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
    
    if (!is.null(offset))
      formMAIN <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                    paste(" + offset(",  offset, ")")) )
    
    model_pois <- glm(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE, family='poisson')
    GoFs_pois <- GoF_stats(model = model_pois)
    
    model_quasipois <- glm(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE, family='quasipoisson')
    GoFs_quasipois <- GoF_stats(model = model_quasipois)
    
    model_negbin <- MASS::glm.nb(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE)
    GoFs_negbin <- GoF_stats(model = model_negbin)
    
    model_zinfl_poisson <- pscl::zeroinfl(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE, dist='poisson')
    GoFs_zinfl_poisson <- GoF_stats(model = model_zinfl_poisson)
    
    model_zinfl_negbin <- pscl::zeroinfl(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE, dist='negbin')
    GoFs_zinfl_negbin <- GoF_stats(model = model_zinfl_negbin)
    
    model_hurdle_poisson <- pscl::hurdle(formMAIN, data=donnes, model=TRUE, x=TRUE, y=TRUE, dist='poisson')
    GoFs_hurdle_poisson <- GoF_stats(model = model_hurdle_poisson)
    
    model_hurdle_negbin <- pscl::hurdle(formMAIN, data=donnes, model=TRUE, x=TRUE, y=TRUE, dist='negbin')
    GoFs_hurdle_negbin <- GoF_stats(model = model_hurdle_negbin)
    
    GoFs <- rbind(GoFs_pois, GoFs_quasipois, GoFs_negbin, GoFs_zinfl_poisson, 
                  GoFs_zinfl_negbin, GoFs_hurdle_poisson, GoFs_hurdle_negbin)
    
    rownames(GoFs) <- c('Poisson','quasi-Poisson','negative binomial',
                        'zero-inflated Poisson','zero-inflated negative binomial','hurdle Poisson','hurdle negative binomial')
    
    if (verbose) {
      message('\nGoodness of Fit for six model types:\n')
      print(round_boc(GoFs, 2), print.gap=4)
      
      # # compare results with vcdExtra::LRstats
      # print(LRstats(model_pois, model_quasipois, model_negbin, 
      # model_zinfl_poisson, model_zinfl_negbin, model_hurdle_poisson, model_hurdle_negbin))
    }
  }
  
  
  # add, if any, factor variables in data to modeldata 
  # (because lm changes names & types and the original variables are needed for PLOTMODEL)
  # factor_variables <- names(modelMAIN$model[sapply(modelMAIN$model, is.factor)])
  factor_variables <- list_xlevels <- names(modelMAIN$xlevels)
  if (!is.null(factor_variables))  modeldata[factor_variables] <- donnes[,factor_variables]
  
  
  output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum,
                 modeldata=modeldata, modelcoefs=modelcoefs,
                 collin_diags=collin_diags, family=family, kind=kind, GoFs=GoFs)
  
  class(output) <- "COUNT_REGRESSION"
  
  return(invisible(output))
  
}


