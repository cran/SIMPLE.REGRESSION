




BF_interps <- function(BF_10 = NULL, BF_01 = NULL, BF_M1 = NULL, BF_M2 = NULL) {
  
  # null & alternative
  if (is.null(BF_10) & !is.null(BF_01))  BF_10 <- 1 / BF_01
  if (is.null(BF_01) & !is.null(BF_10))  BF_01 <- 1 / BF_10
  
  if (!is.null(BF_10) & !is.null(BF_01)) {
    
    message('\nBayes Factor for the null vs. alternative models: ', round(BF_01,2))
    message('\nBayes Factor for the alternative vs. null models: ', round(BF_10,2))
    
    message('\nThe Bayes factor suggests that the observed data are ', round(BF_01,2), ' times more likely')
    message('to occur under the null hypothesis than under the alternative hypothesis.')
    message('\nAlternatively, the observed data are ', round(BF_10,2), ' times more likely to occur')
    message('under the alternative hypothesis than under the null hypothesis.')
    
    # 2016 - Diagnosing the Misuse of the Bayes Factor in Applied Research p 4
    if (BF_01 > BF_10) {
      message('\nSomeone with no prior preference for either hypothesis (i.e., prior odds = 1)')
      message('should now believe that the null model is ', BF_01, ' times more probable')
      message('than the alternative model (i.e., posterior odds = ', BF_01, ' x 1 = ', BF_01, ').')
    }
    
    if (BF_10 > BF_01) {
      message('\nSomeone with no prior preference for either hypothesis (i.e., prior odds = 1)')
      message('should now believe that the alternative model is ', BF_10, ' times more probable')
      message('than the null model (i.e., posterior odds = ', BF_10, ' x 1 = ', BF_10, ').\ n')
    }
    
    if (BF_10 > BF_01) hyp <- 'alternative'
    if (BF_01 > BF_10) hyp <- 'null'
    
    # Jeffreys' scale | Grades or categories of evidence for the Bayes factor
    # Lee M. D. and Wagenmakers, E.-J. (2014) Bayesian cognitive modeling: A 
    # practical course, Cambridge University Press.
    maxBF <- max(c(BF_10, BF_01))
    if (maxBF >= 1  & maxBF < 3)   descr <- 'Anecdotal'
    if (maxBF >= 3  & maxBF < 10)  descr <- 'Moderate'
    if (maxBF >= 10 & maxBF < 30)  descr <- 'Strong'
    if (maxBF >= 30 & maxBF < 100) descr <- 'Very Strong'
    if (maxBF >  100)              descr <- 'Extreme'
    
    message('\nThis level of evidence in favor of the ', hyp, ' hypothesis is') 
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
  
  # prednoms <- names(model.frame(model))[-1]
  
  prednoms <- labels(terms(model))
  
  anova_table <- as.data.frame(anova(model))
  
  lastnom <- prednoms[length(prednoms)]
  
  # for (lupe in 1:(length(prednoms)-1)) {
  # 	
  # 	preds_REORD <- prednoms
  # 	
  # 	preds_REORD[lupe] <- lastnom
  # 	
  # 	preds_REORD[length(prednoms)] <- prednoms[lupe]
  # 		
  # 	form_REORD <- as.formula(paste(DV, paste(preds_REORD, collapse=" + "), sep=" ~ "))
  # 	
  # 	# lm: The terms in the formula will be re-ordered so that main effects come first, 
  # 	# followed by the interactions, all second-order, all third-order and so on: 
  # 	# to avoid this pass a terms object as the formula (see aov and demo(glm.vr) for an example).	
  # 	# model_REORD <- lm(mod_terms, data=data, model=FALSE, x=FALSE, y=FALSE)
  # 	model_REORD <- lm( terms(model), data=data, model=FALSE, x=FALSE, y=FALSE)
  # 	
  # 	anova_table_REORD <- as.data.frame(anova(model_REORD))
  # 
  # 	anova_table[lupe,] <- anova_table_REORD[length(prednoms),]	
  # }	
  
  Eta_Squared <- anova_table$'Sum Sq'[1:length(prednoms)] / 
    (anova_table$'Sum Sq'[1:length(prednoms)] + anova_table$'Sum Sq'[nrow(anova_table)])
  
  anova_table$Eta_Squared <- c(Eta_Squared, NA)
  
  return(invisible(anova_table))
}




diagnostics_plots <- function(modelMAIN, modeldata, plot_diags_nums) {
  
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
    # plot(density(resid(modelMAIN, type='response')))
    
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
      # DV_name <- names(attr(modelMAIN$terms,"dataClasses")[1])
      DV_name <- names(modeldata)[1]
      
      denplotdat_y     <- density(modelMAIN$y)
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
      # same as   plot(modelMAIN, which = 5, sub.caption='')
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
    # cooksd_m0 <- cooks.distance(modelMAIN)  
    # length(cooksd_m0[cooksd_m0 > mean(cooksd_m0) * 2])
    
    # Influence
    # We can also look for influential observations; neither model has any. 
    # which(influence.measures(modelMAIN)$is.inf[,'cook.d'] ) 
    
    if (is.element(13, plot_diags_nums)) {
      # leverage   bar plot-like vertical lines
      # same as   plot(hatvalues(modelMAIN), type = 'h')
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
      # same as    plot(modelMAIN, which = 5, sub.caption='') 
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
      # plot(modelMAIN, which = 6, sub.caption='') 
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

predict_boc <- function(modelMAIN, modeldata, newdata, CI_level = 95, 
                        bootstrap=FALSE, N_sims=100, model_type, family) {
  
  # # doing yhat manually, rather than using built-in predict  prob = factors
  # Xp <- model.matrix(delete.response(terms(modelMAIN)), newdata)
  # b <- coef(modelMAIN)
  # yhat <- c(Xp %*% b)  # c() reshape the single-column matrix to a vector
  
  yhatR <- predict(modelMAIN, newdata=newdata, type="response", se.fit = TRUE)  # type="link"
  
  if (family == 'zinfl_poisson' | family == 'zinfl_negbin') { yhat <- yhatR
  } else { yhat <- yhatR$fit }
  
  
  # CIs
  
  # not doing bootstrapping for MODERATED models because newdata would require values
  # for product terms & lm & glm sometimes alter the var names for interaction terms, too messy
  if (bootstrap & model_type == 'MODERATED') {
    message('\nBootstrapped CIs for a MODERATED regression model was requested but it is not')
    message('available in this version of the package. Regular CI values will be provided instead.\n')
    bootstrap <- FALSE
  }  
  
  
  if (!bootstrap & (family != 'zinfl_poisson' | family != 'zinfl_negbin')) {
    
    # # doing se.fit manually, rather than using built-in predict  prob = factors
    # var.fit <- rowSums((Xp %*% vcov(modelMAIN)) * Xp)  # point-wise variance for predicted mean
    # se.fit <- sqrt(var.fit)
    
    se.fit <- yhatR$se.fit
    
    zforCI <- qnorm((1 + CI_level * .01) / 2)
    
    ci_lb <- yhat - zforCI * se.fit
    
    ci_ub <- yhat + zforCI * se.fit
    
    resmat <- cbind(yhat, ci_lb, ci_ub)
    
    # # not needed when  type="response" on the above predict command
    # if (model_type == 'LOGISTIC')
    # # convert to probabilities
    #   resmat <- 1/(1 + exp(-resmat))
    # 
    # if (model_type == 'COUNT')
    # # exponential function - inverse of log-link
    #   resmat <- exp(resmat)
    
    resmat <- cbind(resmat, se.fit)
  }
  
  if (bootstrap) {
    
    # simulating parameters using the model
    # simulate a range of plausible models based on distribution of parameters from model fit 
    # bootparams <- mvrnorm(n=N_sims, mu=coef(modelMAIN), Sigma=vcov(modelMAIN))
    
    # real bootstrap parameters
    bootparams <- c() 
    
    if (model_type == 'OLS') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- lm(formula(modelMAIN), data=bootsamp)
        
        bootparams <- rbind(bootparams, coef(modelboot))
      }
    }
    
    if (model_type == 'LOGISTIC' | family == 'poisson' |
        family == 'quasipoisson' | family == 'negbin') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modelMAIN$data[sample(rownames(modelMAIN$data), replace=TRUE),]
        
        modelboot <- glm(modelMAIN$formula, family = family, data=bootsamp)
        
        bootparams <- rbind(bootparams, coef(modelboot))
      }
    }
    
    if (family == 'zinfl_poisson') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- pscl::zeroinfl(modelMAIN$formula, dist='poisson', data=bootsamp)
        
        bootparams <- rbind(bootparams, coef(modelboot))
      }
    }
    
    if (family == 'zinfl_negbin') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- pscl::zeroinfl(modelMAIN$formula, dist='negbin', data=bootsamp)
        
        bootparams <- rbind(bootparams, coef(modelboot))
      }
    }
    
    
    ci_lb <- ci_ub <- rep(NA,nrow(newdata))
    
    CI_lo_cut <- (1 - CI_level * .01) / 2
    
    CI_hi_cut <- 1 - (1 - CI_level * .01) / 2
    
    # Xp <- model.matrix(delete.response(terms(modelMAIN)), newdata)   # does not work for factors
    # Xp <- model.frame(delete.response(terms(modelMAIN)), newdata)
    
    Xp <- newdata
    
    # change factors to numeric, use 0 as baseline
    for (lupe in 1:ncol(Xp)) {
      if (is.factor(Xp[,lupe])) {
        Xp[,lupe] <- as.numeric(Xp[,lupe])
        Xp[,lupe] <- 0
      }
    }
    
    # add 1st column of 1s
    Xp <- cbind(1, as.matrix(Xp))
    
    for (lupe in 1:nrow(Xp)) {
      
      if (family == 'zinfl_poisson' | family == 'zinfl_negbin') { simmu <- bootparams %*% cbind(Xp,Xp)[lupe,]
      } else { simmu <- bootparams %*% Xp[lupe,] }
      
      simmu <- sort(simmu)
      
      ci_lb[lupe] <- quantile(simmu, CI_lo_cut)    
      ci_ub[lupe] <- quantile(simmu, CI_hi_cut)
    }
    
    if (model_type == 'LOGISTIC') {
      # convert to probabilities
      ci_lb <- 1/(1 + exp(-ci_lb))   
      ci_ub <- 1/(1 + exp(-ci_ub))   
    }
    
    if (model_type == 'COUNT') {
      # exponential function - inverse of log-link
      ci_lb <- exp(ci_lb)   
      ci_ub <- exp(ci_ub) 
    }
    
    resmat <- cbind(yhat, ci_lb, ci_ub)
  }
  
  
  # 
  # if (!bootstrap) {
  #    
  #   # comparison with built-in predict
  #   yhatR <- predict(modelMAIN, newdata=newdata, type="response", se.fit = TRUE)
  #   print(sum(round(c(resmat[,'yhat']   - yhatR$fit),5)))
  #   print(sum(round(c(resmat[,'se.fit'] - yhatR$se.fit),5)))
  # 
  #   # # comparison with ggeffects
  #   # pgg <- ggeffects::predict_response(modelMAIN, terms = list(AGE = testdata$AGE))
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

