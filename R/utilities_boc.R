


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
    
    new_noms <- c(new_noms, paste('__', paste(comp, base_level, sep='_vs_'), sep=''))
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
    
    # colours for the IV_focal_1_values
    if (length(cols_user) >= length(IV_focal_1_values))  {
      cols <- cols_user[1:length(IV_focal_1_values)]
    } else {cols <- rainbow(length(IV_focal_1_values))}
    
    # if (is.null(ylim)) {
    #   if (CIs) {ylim <- range( pretty(testdata$ci_lb), pretty(testdata$ci_ub))
    #   } else {ylim <- range( pretty(testdata[,DV_predicted]))}
    # }
    # 
    # if (kind != 'LOGISTIC')  ylim[2] <- ylim[2] + (ylim[2] * .05)
    
    if (is.null(IV_focal_2)) {  
      
      don <- t(testdata[DV_predicted]); colnames(don) <- IV_focal_1_values
      
      plot_bar <- barplot(don, beside=TRUE, legend.text=FALSE, col=cols,
                          yaxp = c(min(ylim), max(ylim), 4),
                          ylim=ylim, #yaxt="n", 
                          xpd=FALSE,
                          xlab = IV_focal_1,
                          ylab = ylab,
                          main = title)
      
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
      
      plot_bar <- barplot(plotdon_DV ~ plotdon_IV1 + plotdon_IV2, beside=TRUE,
                          legend.text=FALSE, col=rep(cols, length(IV_focal_1_values)), 
                          yaxp = c(min(ylim), max(ylim), 4),
                          ylim=ylim, #yaxt="n", 
                          xpd=FALSE,
                          xlab = IV_focal_2,
                          ylab = ylab,
                          main = title)
      
      if (CIs) {
        arrows(plot_bar, t(testdata$ci_ub), plot_bar, t(testdata$ci_lb), 
               lwd = 1.0, angle = 90, code = 3, length = 0.05)
      }
      graphics::legend("top", legend = IV_focal_1_values, title = IV_focal_1, 
                       bty="n", lwd=4, col=cols, cex = .80)
    }
  }
}




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
                        bootstrap=FALSE, N_sims=100, kind, family) {
  
  # # computing yhat manually, rather than using built-in predict  NOT USING  prob = factors
  # Xp <- model.matrix(delete.response(terms(modelMAIN)), newdata)
  # b <- coef(modelMAIN)
  # yhat <- c(Xp %*% b)  # c() reshape the single-column matrix to a vector
  
  if (kind == 'ZINFL') {
    yhat_count <- predict(modelMAIN, newdata=newdata, type = "count")
    # zero values are predicted:
    yhat_zero  <- predict(modelMAIN, newdata=newdata, type = "zero", at = 0)  
    # yhat_zero  <- predict(modelMAIN, newdata=newdata, type = "prob")  
  } else if (kind == 'HURDLE') {
      yhat_count <- predict(modelMAIN, newdata=newdata, type = "count")
      # zero values are predicted:
      yhat_zero  <- predict(modelMAIN, newdata=newdata, type = "prob")[,1]
     # yhat_zero_0  <- predict(modelMAIN, newdata=newdata, type = "prob", at = 0:4)
     # print(head(yhat_zero))
     # print(head(yhat_zero_0))
    } else {
      yhatR <- predict(modelMAIN, newdata=newdata, type="response", se.fit = TRUE)
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
      # var.fit <- rowSums((Xp %*% vcov(modelMAIN)) * Xp)  # point-wise variance for predicted mean
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
    # bootparams <- mvrnorm(n=N_sims, mu=coef(modelMAIN), Sigma=vcov(modelMAIN))
    
    # real bootstrap parameters
    predicted_values <- predicted_values_count <- predicted_values_zero <- c()
    
    if (kind == 'OLS') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- lm(formula(modelMAIN), data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }
    
    if (kind == 'LOGISTIC' | kind == 'POISSON') {   
      
      for(i in 1:N_sims) {
        
        bootsamp <- modelMAIN$data[sample(rownames(modelMAIN$data), replace=TRUE),]
        
        modelboot <- glm(modelMAIN$formula, family=family, data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }

    if (kind == 'NEGBIN') {   
      
      for(i in 1:N_sims) {
        
        bootsamp <- modelMAIN$data[sample(rownames(modelMAIN$data), replace=TRUE),]
        
        modelboot <- MASS::glm.nb(modelMAIN$formula, 
                                  family=MASS::negative.binomial(1, link="log"), data=bootsamp)
        
        predicted_values <- cbind(predicted_values, predict(modelboot, newdata=newdata, type="response"))
      }
    }
    
    if (kind == 'ZINFL') {
      
      for(i in 1:N_sims) {
        
        bootsamp <- modeldata[sample(rownames(modeldata), replace=TRUE),]
        
        modelboot <- pscl::zeroinfl(modelMAIN$formula, dist=family, data=bootsamp)
        
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
        
        modelboot <- pscl::hurdle(modelMAIN$formula, dist=family, data=bootsamp)
        
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
  
  