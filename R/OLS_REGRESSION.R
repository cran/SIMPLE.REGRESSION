


OLS_REGRESSION <- function (data, DV, forced=NULL, hierarchical=NULL,
                            COVARS = NULL,
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
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if (!is.null(forced))  {
    donnes <- data[,c(DV,forced)]	
    if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
    formMAIN <- as.formula(paste(DV, paste(forced, collapse=" + "), sep=" ~ "))
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
  if (!is.null(COVARS))        allIVnoms <- c(allIVnoms, COVARS)
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
  
  
  
  # non hierarchical regression (without interactions or blocks)
  if (!is.null(forced))  {
    
    modelMAIN <- lm(formMAIN, data=donnes, model=TRUE, x=TRUE, y=TRUE)
    
    modelMAINsum <- summary(modelMAIN)
    
    RsqMAIN <- modelMAINsum$r.squared
    
    # Anova Table (Type III tests)
    anova_table <- ANOVA_TABLE(data=donnes, model=modelMAIN) 
    
    # creating a new version of the raw data, from the lm output, so that can provide stats for dummy codes
    # modeldata is needed for moderation analyses (next)
    modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
    colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
    
    mainRcoefs <- PARTIAL_COEFS(cormat=cor(modeldata), modelRsq=summary(modelMAIN)$r.squared, verbose=FALSE)
    
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
      
      message('\n\nThe DV is: ', DV)
      message('\nThe IVs are: ', paste(forced, collapse=', ') )
      
      message('\n\nmultiple R = ', round(sqrt(RsqMAIN),3),  
              '   multiple R-squared = ', round(RsqMAIN,3),
              '   adjusted R-squared = ', round(modelMAINsum$adj.r.squared,3))
      
      Fstats <- modelMAINsum$fstatistic
      pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
      
      message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2], '   df_denom = ', Fstats[3], 
              '   p-value = ', round(pvalue,6), '\n')
      
      message('\nAnova Table (Type III tests):\n')
      print(round_boc(anova_table,3), print.gap=4)
      
      modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, confint(modelMAIN))
      message('\n\nModel Coefficients:\n')
      print(round_boc(modelMAINsum$coefficients,3), print.gap=4)
      
      message('\n\nBeta, r, partial correlations, & semi-partial correlations:\n')
      print(round(mainRcoefs,3), print.gap=4)	
      
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
    
    output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum, anova_table=anova_table,
                   mainRcoefs=mainRcoefs, modeldata=modeldata, collin_diags=collin_diags, family='OLS')
    class(output) <- "OLS_REGRESSION"
    
  }
  
  
  
  # hierarchical regression (blocks/steps)
  if (!is.null(hierarchical)) {
    
    for (lupe in 1:length(hierarchical)) {
      
      # keeping info on previous model in lupe for Rsq change stats
      if (lupe > 1) {
        prevModel <- modelMAIN
        prevRsq   <- RsqMAIN
      }
      
      message('\n\n\nStep ', lupe)	
      
      if (lupe==1)  preds <- unlist(hierarchical[1])
      
      if (lupe > 1) preds <- c(preds, unlist(hierarchical[lupe]))
      
      donnesH <- donnes[,c(DV,preds)]
      
      formMAIN <- as.formula(paste(DV, paste(preds, collapse=" + "), sep=" ~ "))
      
      modelMAIN <- lm(formMAIN, data=donnesH, model=TRUE, x=TRUE, y=TRUE)
      
      modelMAINsum <- summary(modelMAIN)
      
      RsqMAIN <- modelMAINsum$r.squared
      
      # Anova Table (Type III tests)
      anova_table <- ANOVA_TABLE(data=donnes, model=modelMAIN) 
      
      # creating a new version of the raw data, from the lm output, so that can provide stats for dummy codes
      # modeldata is needed for moderation analyses (next)
      modeldata <- data.frame(modelMAIN$y, modelMAIN$x[,2:ncol(modelMAIN$x)])
      colnames(modeldata) <- c(DV,colnames(modelMAIN$x)[2:ncol(modelMAIN$x)])
      
      mainRcoefs <- PARTIAL_COEFS(cormat=cor(modeldata), modelRsq=summary(modelMAIN)$r.squared, verbose=FALSE)
      
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
        
        message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2], '   df_denom = ', Fstats[3], 
                '   p-value = ', round(pvalue,6), '\n')
        
        if (lupe > 1) {		
          # current vs previous model comparisons	
          Rsqch <- RsqMAIN - prevRsq
          message('\nCurrent vs previous model comparison:')
          message('\nRsquared change = ', round(Rsqch,3))
          
          fish <- anova(modelMAIN, prevModel)
          message('\nF = ', round(fish$F[2],2), '   df_num = ', 1, '   df_denom = ', fish$Res.Df[1], 
                  '   p-value = ', round(fish$'Pr(>F)'[2],6), '\n')
        }
        
        message('\nAnova Table (Type III tests):\n')
        print(round_boc(anova_table,3), print.gap=4)
        
        message('\nModel Coefficients:\n')
        modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, confint(modelMAIN))
        print(round_boc(modelMAINsum$coefficients,3), print.gap=4)
        
        message('\nBeta, r, partial correlations, & semi-partial correlations:\n')
        print(round(mainRcoefs,3), print.gap=4)	
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
    
    output <- list(modelMAIN=modelMAIN, modelMAINsum=modelMAINsum, anova_table=anova_table, 
                   mainRcoefs=mainRcoefs, modeldata=modeldata, collin_diags=collin_diags, family='OLS')
    class(output) <- "OLS_REGRESSION"
    
  }
  
  
  if (plot_type == 'residuals') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(16,2,3,4))
  
  if (plot_type == 'diagnostics') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(9,12,13,14))
  
  
  return(invisible(output))
  
}


SIMPLE.REGRESSION <- OLS_REGRESSION





