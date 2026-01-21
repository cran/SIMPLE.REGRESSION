

COUNT_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL, formula=NULL,
                              model_type = 'poisson',
                              offset = NULL,
                              CI_level = 95,
                              MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                                  thin = 1, burnin = 1000, 
                                                  HDI_plot_est_type = 'raw'),
                              plot_type = 'residuals',
                              GoF_model_types = TRUE,
                              verbose=TRUE ) {

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
  
  if (!is.null(formula)) {
    formula_vars <- formula_check(formula=formula, data = data)
    DV <- formula_vars$DV
    forced <- formula_vars$forced
    formule <- formula_vars$formule
    hierarchical <- NULL
  }

    if (!is.null(forced))  {
    donnes <- data[,c(DV,forced,offset)]	
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  if (!is.null(hierarchical))  {
    donnes <- data[,c(DV,unlist(hierarchical),offset)]
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  }
  
  # if (is.null(forced) & is.null(hierarchical) ) {
  #   message('\n\nThe data for the analyses were not properly specified.\n')
  #   message('\nThe forced & hierarchical arguments were both NULL (i.e, not provided). Expect errors.\n')
  # }
  
  if (NAflag) cat('\nCases with missing values were found and removed from the data matrix.\n')
  
  
  # gathering all of the predictor variable names
  allIVnoms <- -9999
  if (!is.null(forced))        allIVnoms <- c(allIVnoms, forced)
  if (!is.null(hierarchical))  allIVnoms <- c(allIVnoms, unlist(hierarchical) )
  if (!is.null(offset))        allIVnoms <- c(allIVnoms, offset)
  allIVnoms <- allIVnoms[-1]
  
  donnes <- donnes[c(DV,allIVnoms,offset)]  # a version of donnes that contains only the variables in the analyses
  
  
  # clean up IV factors, contrasts
  IV_type_cleanup(donnes, allIVnoms) 
  
  
  # predictor descriptives
  if (verbose)  descriptives(donnes, allIVnoms)
  
  
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
    # colnames(modelNULLsum$coefs) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    colnames(modelNULLsum$coefs) <- c('B', 'SE', 'z', 'p', 'IRR', 'IRR ci_lb', 'IRR ci_ub')
    
    if (verbose) {
      message('\n\nVariables in the Equation -- NULL model:\n')
      modelNULLsum$coefs <- round_boc(modelNULLsum$coefs, round_non_p = 3, round_p = 6) 
      print(round_boc(modelNULLsum$coefs,3), print.gap=4)
    }
  }
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {    
    
    exp_B <- exp(unlist(modelNULL$coefficients))
    
    exp_B_CIs <- exp(confint.default(modelNULL))
    
    modelNULLsum$coefs$count <- 
      cbind( matrix(modelNULLsum$coefficients$count[1,],nrow=1), exp_B[1], matrix(exp_B_CIs[1,], nrow=1))
    # colnames(modelNULLsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    colnames(modelNULLsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'IRR', 'IRR ci_lb', 'IRR ci_ub')
    
    modelNULLsum$coefs$zero <- 
      cbind( matrix(modelNULLsum$coefficients$zero[1,],nrow=1), exp_B[2], matrix(exp_B_CIs[2,], nrow=1))
    # colnames(modelNULLsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
    colnames(modelNULLsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'IRR', 'IRR ci_lb', 'IRR ci_ub')
    
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
    if (lupe > 1) prevModel <- model
    
    if (verbose) message('\n\n\nBlock ', lupe)	
    
    if (lupe==1)  preds <- unlist(hierarchical[1])
    
    if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
    
    donnesHIER <- donnes[,c(DV,preds,offset)]
    
    # noms_list <- list(DV=DV, IV=preds)

    if (is.null(offset) & is.null(formula))
      formule <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))

    if (!is.null(offset) & is.null(formula))
      formule <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                    paste(" + offset(",  offset, ")")) )

    # if (!is.null(formula)) 
    #   formule <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ ")) 
    
    
        
    if (kind == 'POISSON')
      model <- glm(formule, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, family = family)
    
    if (kind == 'NEGBIN') {
      model <- MASS::glm.nb(formule, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE)
      # model <- glm(formule, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, 
      #                  family = MASS::negative.binomial(1, link="log")) 
    }
    
    if (kind == 'ZINFL')
      model <- pscl::zeroinfl(formule, data = donnesHIER, 
                                  model=TRUE, x=TRUE, y=TRUE, dist = family)
    
    if (kind == 'HURDLE')
      model <- pscl::hurdle(formule, data = donnesHIER, 
                                model=TRUE, x=TRUE, y=TRUE, dist = family)
    
    if (kind == 'ZINFL' | kind == 'HURDLE')
      model$deviance <- -2 * model$loglik
    
    modelsum <- summary(model, correlation=TRUE)
    
    # exponentiated coefficients, & adding them to modelsum$coefs
    
    if (kind == 'POISSON' | kind == 'NEGBIN') {
      exp_B <- exp(model$coefficients)
      exp_B_CIs <- exp(confint.default(model))
      
      modelsum$coefficients <- cbind(modelsum$coefficients, exp_B, exp_B_CIs)
      # colnames(modelsum$coefficients) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
      colnames(modelsum$coefficients) <- c('B', 'SE', 'z', 'p', 'IRR', 'IRR ci_lb', 'IRR ci_ub')
    }
    
    if (kind == 'ZINFL' | kind == 'HURDLE') {
      
      exp_B_count <- exp(model$coefficients$count)
      exp_B_zero  <- exp(model$coefficients$zero)
      
      exp_B_CIs <- exp(confint.default(model))
      
      # dropping the row from the count coefs
      modCcoefs <- modelsum$coefficients$count[row.names(modelsum$coefficients$count) != "Log(theta)",]
      modelsum$coefs$count <- 
        cbind(modCcoefs, exp_B_count, exp_B_CIs[grepl("count", rownames(exp_B_CIs)),])
      # colnames(modelsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
      colnames(modelsum$coefs$count) <- c('B', 'SE', 'z', 'p', 'IRR', 'IRR ci_lb', 'IRR ci_ub')
      
      modelsum$coefs$zero <- 
        cbind(modelsum$coefficients$zero, exp_B_zero, exp_B_CIs[grepl("zero", rownames(exp_B_CIs)),])
      # colnames(modelsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'exp(B)', 'exp(B) ci_lb', 'exp(B) ci_ub')
      colnames(modelsum$coefs$zero) <- c('B', 'SE', 'z', 'p', 'OR', 'OR ci_lb', 'OR ci_ub')
    }
    
    
    # Goodness of Fit
    GoFs <- GoF_stats(model = model)
    
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
    dispersion_ratio_dev <- model$deviance / model$df.residual
    p_deviance <-  1 - pchisq(model$deviance, model$df.residual)
    if (verbose) {
      message('\n\n\nOverdispersion test based on the model deviance:')
      cat('\n    Dispersion Ratio: ', round(dispersion_ratio_dev,3),
          '    Statistic = ', round(model$deviance,3),
          '    p = ', round(p_deviance,5))
    }
    
    # Overdispersion test based on the Pearson Chi-Square
    X2_Pearson <- sum(residuals(model, type = "pearson")^2)
    dispersion_ratio_X2 <- X2_Pearson / model$df.residual
    p_X2 <-  1 - pchisq(X2_Pearson, model$df.residual)
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
    model_dev <- model$deviance 
    # null_dev <- model$null.deviance 
    model_N <- length(model$fitted.values)
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
        print(anova(prevModel, model, test="Chisq"))
      }
      
      if (kind == 'ZINFL' |  kind == 'HURDLE') {
        
        # # LRT to compare these models
        # lmtest::lrtest(model, prevModel)
        # lmtest::waldtest(model, prevModel)
        
        # Compare the current model to a null model without predictors using 
        # chi-squared test on the difference of log likelihoods
        diff_logliks <- as.numeric(2 * (logLik(model) - logLik(prevModel)))
        p <- pchisq(diff_logliks, df = length(preds), lower.tail = FALSE) 
        message('\nComparisons of the current model against the previous model:\n')
        message('\n   Difference in the log likelihoods = ', round(diff_logliks,3), 
                '    df = ', length(preds), '     p = ' , round(p,8), '\n\n')
        
        # Which model should we choose?   https://rpubs.com/lucymark2013/817199
        print(pscl::vuong(model, prevModel))
      }
    }
    
    # Parameter Estimates
    if (verbose) {
      
      if (kind == 'POISSON' | kind == 'NEGBIN') {
        message('\n\nParameter Estimates:\n')
        modelsum$coefs <- round_boc(modelsum$coefficients, round_non_p = 3, round_p = 6) 
        print(round_boc(modelsum$coefs,3), print.gap=4)
        message('\nIRR = incidence rate ratio\n')
      }
      
      if (kind == 'ZINFL' | kind == 'HURDLE') {
        message('\n\nVariables in the Equation -- count portion of the model:\n')
        modelsum$coefs$count <- round_boc(modelsum$coefs$count, round_non_p = 3, round_p = 6) 
        print(round_boc(modelsum$coefs$count,3), print.gap=4)
        
        message('\n\nVariables in the Equation -- zero portion of the model:\n')
        modelsum$coefs$zero <- round_boc(modelsum$coefs$zero, round_non_p = 3, round_p = 6) 
        print(round_boc(modelsum$coefs$zero,3), print.gap=4)
        
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
    
    MCMC_outp <- list(Bayes_HDIs = NULL, chains = NULL, autocorrels = NULL, SDs = NULL)
    
    if (MCMC_options$MCMC & (kind == 'POISSON' | kind == 'NEGBIN')) {

      if (family == 'poisson') famille <- family
      
      if (family == 'quasipoisson') {
        message("\n\nFamily = 'quasipoisson' analyses are currently not possible for")
        message("the MCMC analyses. family = 'poisson' will therefore be used instead.\n")
        famille <- 'poisson'
      }
      
      if (kind == 'NEGBIN') famille <- "neg_binomial_2"
        
      # MCMC_mod <- rstanarm::stan_glm(formule, data = donnesHIER, family = famille,
      #                                refresh = 0, algorithm="sampling", 
      #                                iter = MCMC_options$Nsamples)
      
      MCMC_outp <- Bayes_HDIs_reg_preds(formule, data = donnesHIER, CI_level,  
                                        MCMC_options, #noms_list, 
                                        model_type = 'COUNT', famille = famille)
      
      if (verbose) {
        message('\n\nBayesian MCMC chains:')
        message('\n    Number iterations (samples): ', MCMC_options$Nsamples)
        message('\n    Thinning interval:  ', MCMC_options$thin)
        message('\n    Burn-in period:  ', MCMC_options$burnin)
        message('\n    Final chain length (# of samples): ', nrow(MCMC_outp$chains))
        message('\n    The priors were the default rstanarm::stan_glm function priors')
        
        message('\n\nAutocorrelations for the MCMC chains:\n')
        print(round(MCMC_outp$autocorrels,3), print.gap=4)
        
        message('\n\nBayesian HDIs for the raw coefficients & the exponentiated raw coefficients:\n')
        print(round_boc(MCMC_outp$Bayes_HDIs, round_non_p = 2, round_p = 5), print.gap=4)
      }
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
        
        if (kind == 'POISSON')
          M1 <- glm(formM1, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (kind == 'NEGBIN') {
          # M1 <- MASS::glm.nb(formM1, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE)
          M1 <- glm(formM1, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log")) 
        }
        
        if (kind == 'ZINFL')
          M1 <- pscl::zeroinfl(formM1, data = donnesHIER, 
                               model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (kind == 'HURDLE')
          M1 <- pscl::hurdle(formM1, data = donnesHIER, 
                             model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (is.null(offset))
          formM2 <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
        
        if (!is.null(offset))
          formM2 <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                      paste(" + offset(",  offset, ")")) )
        
        if (kind == 'POISSON')
          M2 <- glm(formM2, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, family = family)
        
        if (kind == 'NEGBIN') {
          # M2 <- MASS::glm.nb(formM2, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE)
          M2 <- glm(formM2, data = donnesHIER, model=TRUE, x=TRUE, y=TRUE, 
                    family = MASS::negative.binomial(1, link="log"))
        }
        
        if (kind == 'ZINFL')
          M2 <- pscl::zeroinfl(formM2, data = donnesHIER, 
                               model=TRUE, x=TRUE, y=TRUE, dist = family)
        
        if (kind == 'HURDLE')
          M2 <- pscl::hurdle(formM2, data = donnesHIER, 
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
      print(suppressWarnings(anova(model, test = 'Rao')))
      
      message('\n\nCorrelations of Parameter Estimates:\n')
      print(round_boc(modelsum$correlation,3), print.gap=4)
    }
  }
  
  
  # modeldata
  if (kind == 'POISSON')  modeldata <- model$data
  
  if (kind == 'NEGBIN') {
    
    if (is.null(offset)) {
      
      modeldata <- data.frame(model$y, model$x[,2:ncol(model$x)])
      
      colnames(modeldata) <- c(DV,colnames(model$x)[2:ncol(model$x)])
    }
    
    if (!is.null(offset)) {
      
      modeldata <- data.frame(model$y, model$x[,2:ncol(model$x)], model$offset)
      
      colnames(modeldata) <- c(DV,colnames(model$x)[2:ncol(model$x)], offset)
    }
  }
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    # modeldata <- data.frame(model$y, model$x$zero[,2:ncol(model$x$zero)])
    # colnames(modeldata) <- c(DV,colnames(model$x$zero)[2:ncol(model$x$zero)])
    modeldata <- model$model
    
    if (!is.null(offset))  colnames(modeldata)[ncol(modeldata)] <- offset
  }
  
  modeldata$predicted <- model$fitted.values
  
  
  # modelcoefs <- modelsum$coefs
  modelcoefs <- modelsum$coefficients
  
  
  # variable correlations - using model data with dummy coded factor variables (all numeric)
  if (kind == 'POISSON' | kind == 'NEGBIN') {
    modeldataN <- data.frame(model$y, model$x[,2:ncol(model$x)])
    colnames(modeldataN) <- c(DV,colnames(model$x)[2:ncol(model$x)])
  }
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    modeldataN <- data.frame(model$y, model$x$zero[,2:ncol(model$x$zero)])
    colnames(modeldataN) <- c(DV,colnames(model$x$zero)[2:ncol(model$x$zero)])
  }
  modeldatacorrels <- cor(modeldataN)
  
  
  # casewise diagnostics
  modeldata$residuals_raw <-          resid(model, type='response')
  modeldata$residuals_pearson <-      resid(model, type='pearson')
  
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
    modeldata$residuals_pearson_std <-  rstandard(model, type='pearson') 
    modeldata$residuals_deviance <-     resid(model, type='deviance')
    modeldata$residuals_deviance_std <- rstandard(model, type='deviance')
    modeldata$residuals_quantile <-     quantile_residuals(model)
    modeldata$cooks_distance <-         cooks.distance(model)
    modeldata$dfbeta <-                 dfbeta(model)
    modeldata$dffit <-                  dffits(model)
    modeldata$leverage <-               hatvalues(model)
    modeldata$covariance_ratios <-      covratio(model)
  }
  
  collin_diags <- Collinearity(model.matrix(model), verbose=FALSE)
  
  if (verbose) {			
    message('\n\nPredictor variable correlation matrix:\n')
    print(round(modeldatacorrels,3), print.gap=4)
    
    message('\n\nCollinearity Diagnostics:')
    VIFtol_messages(collin_diags$VIFtol[,2])
    CondInd_messages(collin_diags$CondInd)
  }     

  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    
    diagnostics_plots(model=model, modeldata=modeldata, plot_diags_nums=c(5, 9))
    
  } else {
    
    if (plot_type == 'residuals') 
      
      diagnostics_plots(model=model, modeldata=modeldata, plot_diags_nums=c(8, 9, 10, 11))
    
    if (plot_type == 'diagnostics')
      
      diagnostics_plots(model=model, modeldata=modeldata, plot_diags_nums=c(12, 13, 14, 15))
    
    if (plot_type == 'Bayes_HDI' & MCMC_options$MCMC & !is.null(MCMC_outp$chain)) 
      
      Bayes_HDI_plot(chain_dat = MCMC_outp$chains, 
                     Bayes_HDIs = MCMC_outp$Bayes_HDIs,
                     CI_level = CI_level, # SDs = xx$SDs[c(preds)],
                     HDI_plot_est_type = MCMC_options$HDI_plot_est_type)
  }
  
  
  
  if (GoF_model_types) {
    # Goodness of Fit for all possible model types
    
    preds <- unlist(hierarchical[1])
    
    donnesHIER <- donnes[,c(DV,preds,offset)]
    
    if (is.null(offset))
      formule <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
    
    if (!is.null(offset))
      formule <- as.formula( paste(paste(DV, paste(preds, collapse=" + "), sep=" ~ "), 
                                    paste(" + offset(",  offset, ")")) )
    
    model_pois <- glm(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, family='poisson')
    GoFs_pois <- GoF_stats(model = model_pois)
    
    model_quasipois <- glm(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, family='quasipoisson')
    GoFs_quasipois <- GoF_stats(model = model_quasipois)
    
    model_negbin <- MASS::glm.nb(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE)
    GoFs_negbin <- GoF_stats(model = model_negbin)
    
    model_zinfl_poisson <- pscl::zeroinfl(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, dist='poisson')
    GoFs_zinfl_poisson <- GoF_stats(model = model_zinfl_poisson)
    
    model_zinfl_negbin <- pscl::zeroinfl(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, dist='negbin')
    GoFs_zinfl_negbin <- GoF_stats(model = model_zinfl_negbin)
    
    model_hurdle_poisson <- pscl::hurdle(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, dist='poisson')
    GoFs_hurdle_poisson <- GoF_stats(model = model_hurdle_poisson)
    
    model_hurdle_negbin <- pscl::hurdle(formule, data=donnesHIER, model=TRUE, x=TRUE, y=TRUE, dist='negbin')
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
  
  
  # add, if any, factor variables in data to modeldata, for PLOTMODEL)
  # factor_variables <- names(model$model[sapply(model$model, is.factor)])
  factor_variables <- list_xlevels <- names(model$xlevels)
  if (!is.null(factor_variables))  modeldata[factor_variables] <- donnes[,factor_variables]
  
  
  output <- list(model=model, modelsum=modelsum,
                 modeldata=modeldata, modelcoefs=modelcoefs,
                 collin_diags=collin_diags, family=family, kind=kind, GoFs=GoFs,
                 chain_dat = MCMC_outp$chains,
                 Bayes_HDIs = MCMC_outp$Bayes_HDIs) # noms_list = noms_list)
  
  class(output) <- "COUNT_REGRESSION"
  
  return(invisible(output))
  
}


