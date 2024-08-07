

MODERATED_REGRESSION <- function (data, DV, IV, MOD,
                                  IV_type = 'numeric', IV_range='tumble',
                                  MOD_type = 'numeric', MOD_levels='quantiles', MOD_range=NULL,
                                  quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                                  COVARS = NULL,
                                  center = TRUE,  
                                  plot_type = 'residuals', plot_title=NULL, DV_range=NULL,
                                  Xaxis_label=NULL, Yaxis_label=NULL, legend_label=NULL,
                                  JN_type = 'Huitema',
                                  verbose=TRUE ) {
  
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if ( is.null(IV) | is.null(MOD) ) {
    message('\n\nThe data for the analyses were not properly specified.\n')
    message('\nThe IV or MOD arguments were NULL. Expect errors.\n')
  }
  
  if (is.numeric(IV_range)) {
    IV_range_user <- IV_range
    IV_range <- 'numeric'
  }
  
  if (is.factor(data[,MOD]))  MOD_type <- 'factor'
  
  if (MOD_type == 'factor' & !is.factor(data[,MOD])) data[,MOD] <- factor(data[,MOD])
  
  donnes <- data[,c(DV,IV,MOD,COVARS)]
  
  if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  
  formMAIN <- as.formula(paste(DV, paste(c(IV,MOD,COVARS), collapse=" + "), sep=" ~ "))
  
  
  if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
  
  
  if (length(MOD_levels) > 1) {
    
    mod_values <- MOD_levels
    
    MOD_levels <- 'user'
  }
  
  
  # IV_type
  if (IV_type == 'numeric' & is.numeric(donnes[,IV])) IV_type <- 'numeric'
  
  # making IV a factor if it is non numeric & has 2 levels
  if (!is.numeric(donnes[,IV])) {
    if (is.factor(donnes[,IV])) {
      if (length(levels(donnes[,IV])) == 2) 
        message('\n\nThe IV in data is a factor with two levels.')
      if (length(levels(donnes[,IV])) > 2) {
        message('\n\nThe IV in data is a factor with > 2 levels.')
        message('\n\nPlots cannot be produced. The analyses cannot be performed.')
      }
    }
    if (!is.factor(donnes[,IV])) {
      donnes[,IV] <- factor(donnes[,IV])
      if (length(levels(donnes[,IV])) == 2) {
        message('\n\nThe IV in data is non numeric and has been converted to a')
        message('\nfactor with two levels.')
      }
      if (length(levels(donnes[,IV])) > 2) {
        message('\n\nThe IV in data is non numeric, it has been converted to a')
        message('\nfactor, but there are > 2 levels.')
        message('\n\n Plots cannot be produced. The analyses cannot be performed.')
      }		
    }
    IV_type <- 'factor'
  }
  
  
  # centering, if requested
  if (center) {
    if (is.numeric(donnes[,IV]))   donnes[,IV]  <- donnes[,IV]  - mean(donnes[,IV])	
    if (is.numeric(donnes[,MOD]))  donnes[,MOD] <- donnes[,MOD] - mean(donnes[,MOD])	
  }
  
  
  # gathering all of the predictor variable names
  allIVnoms <- c(IV, MOD)
  if (!is.null(COVARS))  allIVnoms <- c(allIVnoms, COVARS)
  
  # a version of donnes that contains only the variables in the analyses
  donnesRED <- donnes[c(DV,allIVnoms)]  
  
  
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
  

  
  modelMAIN <- lm(formMAIN, data=donnes, model=TRUE, x=TRUE, y=TRUE)
  
  modelMAINsum <- summary(modelMAIN)
  
  RsqMAIN <- modelMAINsum$r.squared
  
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
      message('\nThe IV is: ', IV)
      message('\nThe MOD is: ', MOD)
      message('\n\nMultiple regression statistics for when the product term is not in the model:')

    message('\n\nmultiple R = ', round(sqrt(RsqMAIN),3),  
            '   multiple R-squared = ', round(RsqMAIN,3),
            '   adjusted R-squared = ', round(modelMAINsum$adj.r.squared,3))
    
    Fstats <- modelMAINsum$fstatistic
    pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
    
    message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2], '   df_denom = ', Fstats[3], 
            '   p-value = ', round(pvalue,6), '\n')
    
    modelMAINsum$coefficients <- cbind(modelMAINsum$coefficients, confint(modelMAIN))
    message('\nModel Coefficients:\n')
    print(round_boc(modelMAINsum$coefficients,3), print.gap=4)
    
    message('\nBeta, r, partial correlations, & semi-partial correlations:\n')
    print(round(mainRcoefs,3), print.gap=4)	
    
    message('\n\nCorrelation matrix:\n')
    # if (!is.null(forced)) print(round(cor(modeldata[,c(DV,forced)]),3), print.gap=4)		
    # if (modreg)           print(round(cor(modeldata[,c(IV,MOD,COVARS)]),3), print.gap=4)
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
  
  
  
  
  # get the names of the MOD variable/its dummy codes
  noms <- colnames(modeldata)
  
  # names of the regression diagnostic variables -- for exclusion from MODnew
  nomsdiags <- c("residuals_standardized", "residuals_studentized", "cooks_distance",        
                 "dfbeta", "dffit", "leverage", "covariance_ratios", "predicted", "residuals")
  
  MODnew <- setdiff(noms,c(DV,IV,COVARS, nomsdiags))
  
  
  # create PROD(s)
  prods <- data.frame(modeldata[,IV] * modeldata[,MODnew])
  colnames(prods) <- paste(paste(IV, ":", sep = ""), MODnew, sep = "")
  head(prods)
  
  PROD <- colnames(prods)
  don5 <- data.frame(cbind(modeldata, prods), check.names=FALSE)
  
  termsPOSSIBLE <- c(IV,MODnew,PROD,COVARS)
  
  formXn <- as.formula(paste(DV, paste(termsPOSSIBLE, collapse=" + "), sep=" ~ "))
  
  modelXN <- lm(formXn, data=don5, model=TRUE, x=TRUE, y=TRUE)
  #print(modelXN); print(summary(modelXN)); print(vcov(modelXN))
  
  modelXNsum <- summary(modelXN)
  
  RsqXN <- modelXNsum$r.squared
  
  # XN vs MAIN model comparisons	
  RsqchXn <- RsqXN - RsqMAIN	
  fsquaredXN <- (RsqXN - RsqMAIN) / (1 - RsqXN)	
  
  xnRcoefs <- PARTIAL_COEFS(cormat=cor(don5[,c(IV,MODnew,PROD,COVARS)]), modelRsq=summary(modelXN)$r.squared, verbose=FALSE)
  
  modeldata$predicted <- modelXN$fitted.values
  
  
  # casewise diagnostics
  modeldata$residuals <- resid(modelXN)
  modeldata$residuals_standardized <- rstandard(modelXN)
  modeldata$residuals_studentized <- rstudent(modelXN)
  modeldata$cooks_distance <- cooks.distance(modelXN)
  modeldata$dfbeta <- dfbeta(modelXN)
  modeldata$dffit <- dffits(modelXN)
  modeldata$leverage <- hatvalues(modelXN)
  modeldata$covariance_ratios <- covratio(modelXN)
  
  
  if (verbose) {	
    
    message('\n\nModerated regression:')
    
    if (MOD_type == 'numeric' & MOD_levels == 'user')
      message('\nThe specification for MOD_levels is: ', paste(mod_values, collapse=', '), '\n')
    
    if (MOD_levels == 'quantiles') 
      message('\nThe moderator quantiles that will be used: ', paste(quantiles_MOD, collapse=', ') )
    
    message('\n\nmultiple R = ', round(sqrt(RsqXN),3), 
            '   multiple R-squared = ', round(RsqXN,3),
            '   adjusted R-squared = ', round(modelXNsum$adj.r.squared,3))
    
    Fstats <- modelXNsum$fstatistic
    pvalue <- pf(Fstats[1], Fstats[2], Fstats[3], lower.tail=FALSE)
    
    message('\nF = ', round(Fstats[1],2), '   df_num = ', Fstats[2], '   df_denom = ', Fstats[3], 
            '   p-value = ', round(pvalue,6))
    
    message('\n\nModel comparisons:')
    message('\nRsquared change = ', round(RsqchXn,3))
    message('\nf-squared change = ', round(fsquaredXN,3))
    
    fish <- anova(modelXN, modelMAIN)
    message('\nF = ', round(fish$F[2],2), '   df_num = ', 1, '   df_denom = ', fish$Res.Df[1], 
            '   p-value = ', round(fish$'Pr(>F)'[2],6), '\n')
    
    modelXNsum$coefficients <- cbind(modelXNsum$coefficients, confint(modelXN))
    message('\nModel Coefficients:\n')
    print(round(modelXNsum$coefficients,3), print.gap=4)
    
    message('\nBeta, r, partial correlations, & semi-partial correlations:\n')
    print(round(xnRcoefs,3), print.gap=4)
  }	
  
  
  
  # the values/levels of MOD, if it is continuous
  if (MOD_type != 'factor') {  
    
    if (MOD_levels=='quantiles')  mod_values <- quantile(donnes[,MOD], probs=quantiles_MOD)
    
    if (MOD_levels=='AikenWest') {
      MODmn <- sapply(donnes[MOD], mean, na.rm = TRUE)
      MODsd <- sapply(donnes[MOD], sd, na.rm = TRUE)
      MODlo <- MODmn - MODsd
      MODhi <- MODmn + MODsd
      mod_values <- c(MODlo, MODmn, MODhi)
    }	
    
    # if (length(MOD_levels) > 1) mod_values <- MOD_levels
    
    # mod_values <- round(mod_values,3)
    
    # make sure there are no duplicates, due to rounding, in mod_values
    for (lupe in 1:10) {
      if (!any(duplicated(round(mod_values,lupe)))) { 
        mod_values <- round(mod_values,lupe)
        break
      }
    }
  }
  
  # the values/levels of MOD, if it is categorical
  if (MOD_type == 'factor')  mod_values <- levels(donnes[,MOD])
  
  
  
  
  
  ###########################  simple slopes  ###########################
  
  coefs <- modelXN$coefficients
  Sb <- vcov(modelXN)[2:length(coefs),2:length(coefs)]
  if (MOD_type != 'factor') {
    slopes   <- coefs[IV] + coefs[PROD] * mod_values
    intercepts  <- coefs[MOD] * mod_values + coefs['(Intercept)']
    SEslopes <- sqrt( Sb[IV,IV] + 2*mod_values*Sb[IV,PROD] +  mod_values**2 * Sb[PROD,PROD])
  } else { 
    slopes   <- c(coefs[IV], (coefs[IV] + coefs[PROD]))
    intercepts  <- c(coefs[1],  (coefs[1]  + coefs[MODnew]))
    diagSB   <- diag(Sb)
    SEslopes <- c(sqrt(Sb[1,1]), sqrt( abs(Sb[1,1] - diagSB[PROD])))
  }
  tslopes <- slopes / diag(SEslopes)
  tslopes <- slopes / SEslopes
  N <- nrow(donnes)
  df <- modelXN$df.residual
  k <- length(coefs) - 1                    # number of predictors    # FIX?
  dfs <-  N - k -1  
  pslopes <- (1 - pt(abs(tslopes),dfs)) * 2
  # CIs - 2003 Cohen Aiken West p 278
  tabledT <- qt(.975, dfs)
  me <- tabledT * SEslopes
  confidLo <- slopes - me
  confidHi <- slopes + me
  
  # simslop <- data.frame(MODlevel=mod_values, Intercept=intercepts, Slope=slopes, SEslopes=SEslopes, 
  # t=tslopes, p=pslopes, Slope_CI_lo=confidLo, Slope_CI_hi=confidHi)
  simslop <- data.frame(Intercept=intercepts, b=slopes, SE_b=SEslopes, 
                        t=tslopes, p=pslopes, b_CI_lo=confidLo, b_CI_hi=confidHi)
  rownames(simslop) <- mod_values
  
  if (verbose) {
    message('\n\nSimple slopes for the levels of the moderator:\n')
    print(round_boc(simslop), print.gap=4)   # , row.names = FALSE)
  }
  
  # standardized slopes
  if (MOD_type == 'factor') {	
    grp_correls <- sapply(
      split(data.frame(donnes[,DV], donnes[,IV]), donnes[,MOD]),     # FIX?
      function(x) cor(x[[1]],x[[2]]) )
    grp_DV_SDs <- tapply(donnes[,DV], donnes[,MOD], sd)
    grp_IV_SDs <- tapply(donnes[,IV], donnes[,MOD], sd)
    # grp_DV_MNs <- tapply(donnes[,DV], donnes[,MOD], mean)   # for JN
    grp_IV_MNs <- tapply(donnes[,IV], donnes[,MOD], mean)   # for JN
    zslopes <- slopes   * (grp_IV_SDs / grp_DV_SDs)
    zSE     <- SEslopes * (grp_IV_SDs / grp_DV_SDs)
    me <- tabledT * zSE
    confidLo <- zslopes - me
    confidHi <- zslopes + me
    simslopZ <- data.frame(beta=zslopes, SE_beta=zSE, 
                           beta_CI_lo=confidLo, beta_CI_hi=confidHi, r=grp_correls)
    rownames(simslopZ) <- mod_values
  } else {
    zslopes <- slopes   * (sd(modeldata[,IV]) / sd(modeldata[,DV]) )
    zSE     <- SEslopes * (sd(modeldata[,IV]) / sd(modeldata[,DV]) )	
    # compute zslopes  = slopes &*  (sd(1,1)/sd(1,4)).
    # compute zSE = SEslopes &* (sd(1,1)/sd(1,4))  .	
    me <- tabledT * zSE
    confidLo <- zslopes - me
    confidHi <- zslopes + me	
    reffsize <- sqrt( tslopes**2 / (tslopes**2 + dfs) )  # EffectSizeConversion.pdf	
    simslopZ <- data.frame(beta=zslopes, SE_beta=zSE, 
                           beta_CI_lo=confidLo, beta_CI_hi=confidHi, r=reffsize)
    rownames(simslopZ) <- mod_values
  }
  
  if (verbose) {
    message('\nStandardized simple slopes & r for the levels of the moderator:\n')
    print(round_boc(simslopZ), print.gap=4)   # , row.names = FALSE)
  }
  
  
  
  
  ###########################  Johnson-Neyman regions of significance  ###########################
  
  IV_min_JN   <- min(modelXN$model[,IV])
  IV_max_JN   <- max(modelXN$model[,IV])
  
  # # IV min & max 
  # if (is.null(IV_range)) {
  # IV_min_JN <- min(modelXN$x[,IV])
  # IV_max_JN <- max(modelXN$x[,IV])
  # } else {
  # IV_min_JN <- IV_range[1]
  # IV_max_JN <- IV_range[2]
  # }		
  
  # modelXN  <- model$modelXN
  
  JN.data <- NULL
  
  ros <- NULL
  
  if (MOD_type == 'numeric') {
    
    # MOD values 
    if (is.null(MOD_range)) {
      MOD_min <- min(modelXN$x[,MOD])
      MOD_max <- max(modelXN$x[,MOD])
    } else {
      MOD_min <- MOD_range[1]
      MOD_max <- MOD_range[2]
    }
    byint <- abs(MOD_min - MOD_max) / 100
    MODvalues <- seq(MOD_min, MOD_max, by = byint)
    
    
    B <- modelXN$coefficients[c('(Intercept)',IV,MOD,PROD)]
    S <- stats::vcov(modelXN)[c('(Intercept)',IV,MOD,PROD), c('(Intercept)',IV,MOD,PROD)]
    
    tcrit <- qt(0.975, modelXN$df.residual)     
    
    # construct the quadratic equation
    a <- tcrit^2 * S[PROD,PROD] - B[PROD]^2
    b <- 2 * (tcrit^2 * S[IV,PROD] - B[IV]*B[PROD])
    c <- tcrit^2 * S[IV,IV] - B[IV]^2
    
    
    # when JN values can be computed
    if (b^2-4*a*c >= 0) {
      JNa <- (-b - sqrt(b^2-4*a*c)) / (2*a)
      JNb <- (-b + sqrt(b^2-4*a*c)) / (2*a)
      ros <- sort(c(JNa,JNb))
      
      # adding the ros values to MODvalues, but only if they are within the MOD min/max range
      MODvalues <- sort(c(MODvalues, ros[ros >= MOD_min & ros <= MOD_max ])	)
      
      # data for plot
      SimpleSlope <- B[IV]+(B[PROD]*MODvalues)
      StdError <- sqrt((S[IV,IV])+(MODvalues^2*S[PROD,PROD])+(2*MODvalues*S[IV,PROD]) )   
      CI.L <- SimpleSlope - tcrit*StdError  
      CI.U <- SimpleSlope + tcrit*StdError  
      JN.data <- data.frame(SimpleSlope, CI.L, CI.U, MODvalues, StdError)
      
      if (verbose) {
        
        message('\n\nJohnson-Neyman regions of significance interval:')  
        message('\n   The slopes for ',IV,' predicting ',DV,' are significant when ',MOD)
        message('\n','   values are outside of (lower and higher than) this range: ',round(ros[1],2),' to ',round(ros[2],2))
        
        # slope & CI values at the ros points
        if (ros[1] >= MOD_min) {
          rownum <- which(JN.data$MODvalues == ros[1]) # the row with the low ros MOD value info
          message('\n   When ',MOD,' is ',round(ros[1],2),
                  ':   slope = ',round(JN.data$SimpleSlope[rownum],2),
                  '  CI_lb = ',round(JN.data$CI.L[rownum],2),'  CI_ub = ',round(JN.data$CI.U[rownum],2),
                  '  SE = ',round(JN.data$StdError[rownum],2))
        }
        if (ros[2] <= MOD_max) {
          rownum <- which(JN.data$MODvalues == ros[2]) # the row with the high ros MOD value info
          message('\n   When ',MOD,' is ',round(ros[2],2),
                  ':    slope = ',round(JN.data$SimpleSlope[rownum],2),
                  '    CI_lb = ',round(JN.data$CI.L[rownum],2),'   CI_ub = ',round(JN.data$CI.U[rownum],2),
                  '    SE = ',round(JN.data$StdError[rownum],2))
        }
        message('\n   The ',IV,' values range from ',round(IV_min_JN,2),' to ',round(IV_max_JN,2), '\n\n')
      }
    }	
    
    if (b^2-4*a*c < 0) {
      message('\n\nJohnson-Neyman regions of significance cannot be computed because the')  
      message('procedure would involve finding the square root of a negative number.\n')
    }		
  }
  
  
  
  if (MOD_type == 'factor') {
    
    N <- length(modelXN$y)
    Ngroups <- length(MODnew) + 1  # OK?	 	
    Ncombins <- choose(Ngroups,2)  # the # of possible 2-group comparisons
    grpNs <- grp_IV_MNs <- grp_DV_MNs <- sscp11 <- sscp12 <- sscp22 <- ssreg <- 
      intercepts <- slopes <- ssresid <- rep(-9999, Ngroups)
    JNlohigrps <- c(NA,NA)
    compnoms <- NA
    # Xhi <- Xlo <- NA  # matrix(NA, choose(Ngroups,2), 2) 
    alpha <- .05  # FIX
    
    MOD_grps <- levels(donnes[,MOD])
    
    for (lupe in 1:Ngroups) {
      
      don <- (subset(donnes, donnes[MOD] == MOD_grps[lupe], select=c(IV,DV)))
      
      don2 <- data.frame(IV=unlist(don[IV]), DV = unlist(don[DV]))
      colnames(don2) <- c(IV,DV)
      
      formSIMPLE <- as.formula(paste(DV, IV, sep=" ~ "))
      
      donmod <- lm(formSIMPLE, data=don2)
      
      intercepts[lupe] <- donmod$coefficients['(Intercept)']
      slopes[lupe]  <- donmod$coefficients[IV]
      
      grpNs[lupe] <- nrow(don)
      
      grp_IV_MNs[lupe] <- mean(don[,1])
      grp_DV_MNs[lupe] <- mean(don[,2])
      
      sscp <- cov(don) * (grpNs[lupe] - 1)
      
      sscp11[lupe] <- sscp[1,1]
      sscp12[lupe] <- sscp[1,2]
      sscp22[lupe] <- sscp[2,2]
      
      ssreg[lupe] <- (sscp12[lupe]**2) / sscp11[lupe]
      
      sumdonmod <- summary(donmod)
      ssresid[lupe] <- sum(sumdonmod$residuals**2)
    }
    
    # Huitema 1981 p 293,   also Huitema 2011 p 277
    bigC <- Ncombins
    bigJ <- Ngroups
    bigCprime <- (bigJ * (bigJ - 1)) / 2  # same as Ncombins
    
    if (JN_type == 'Huitema') {
      
      if (Ngroups == 2) {			
        Fcrit <- qf(p=alpha, df1 =  1,       df2 = (N - bigJ * (bigC + 1)), lower.tail=FALSE) 
      }
      if (Ngroups > 2) {
        alpha = 0.05 / bigCprime
        # If more than two groups are involved, the Bonferroni F is    (Huitema 1981 p 293)
        Fcrit <- qf(p=alpha, df1 = bigC + 1, df2 = (N - bigJ * (bigC + 1)), lower.tail=FALSE) 
      }	
    }		
    
    if (JN_type == 'Pedhazur') 
      Fcrit <- qf(p=alpha, df1 =  2, df2 = N - 4, lower.tail=FALSE) 
    
    
    # comparing groups
    for (lupe1 in 1:(Ngroups-1)) {
      for (lupe2 in (lupe1+1):(Ngroups)) {
        
        ssresd <- sum(ssresid) # works for Pedhazur
        
        if (JN_type == 'Huitema') {  # Huitema 2011 p 253
          
          A <- (-Fcrit / (N - 4)) * ssresd * 
            (1/sscp11[lupe1] + 1/sscp11[lupe2]) + (slopes[lupe1] - slopes[lupe2])**2
          
          B <- ( Fcrit / (N - 4)) * ssresd * (grp_IV_MNs[lupe1]/sscp11[lupe1] + 
                                                grp_IV_MNs[lupe2]/sscp11[lupe2]) + 
            (intercepts[lupe1] - intercepts[lupe2]) * (slopes[lupe1] - slopes[lupe2])
          
          C <- (-Fcrit / (N - 4)) * ssresd * ( (N/(grpNs[lupe1]*grpNs[lupe2])) + 
                                                 (grp_IV_MNs[lupe1]**2/sscp11[lupe1]) + 
                                                 (grp_IV_MNs[lupe2]**2/sscp11[lupe2]) ) + (intercepts[lupe1] - intercepts[lupe2])**2
        }
        
        
        if (JN_type == 'Pedhazur') {   # Pedhazur 1997 p 593 / Potoff 1964
          A <- (-2 * Fcrit / (N - 4)) * ssresd * (1/sscp11[lupe1] + 1/sscp11[lupe2]) + 
            (slopes[lupe1] - slopes[lupe2])**2
          
          B <- ( 2 * Fcrit / (N - 4)) * ssresd * 
            (grp_IV_MNs[lupe1]/sscp11[lupe1] + grp_IV_MNs[lupe2]/sscp11[lupe2]) + 
            (intercepts[lupe1] - intercepts[lupe2]) * (slopes[lupe1] - slopes[lupe2])
          
          C <- (-2 * Fcrit / (N - 4)) * ssresd * 
            ( (N/(grpNs[lupe1]*grpNs[lupe2])) + (grp_IV_MNs[lupe1]**2/sscp11[lupe1]) + 
                (grp_IV_MNs[lupe2]**2/sscp11[lupe2]) ) + (intercepts[lupe1] - intercepts[lupe2])**2				     		
        }		
        
        # XL1 <- (B*-1 - sqrt(B**2 - A * C)) / A
        
        # XL2 <- (B*-1 + sqrt(B**2 - A * C)) / A
        
        # lohi <- rbind(lohi, c(XL1, XL2))
        
        if ((B**2 - A*C) > 0) {
          hi <- (-B + (sqrt(B**2 - A*C)) ) / A
          lo <- (-B - (sqrt(B**2 - A*C)) ) / A
          
          if (hi > lo)  JNlohigrps <- rbind(JNlohigrps, c(lo, hi))
          
          # if (hi > lo)  Xhi <- rbind(Xhi, hi)
          # if (lo < hi)  Xlo <- rbind(Xlo, lo)
        } else {
          JNlohigrps <- rbind(JNlohigrps, c(NA,NA))
        }	
        
        
        compnoms <- rbind(compnoms,  paste('Groups',lupe1,'&',lupe2, sep=' ') )
      }
    }
    
    JNlohigrps <- JNlohigrps[-1, ,drop = FALSE]
    # JNlohigrps <- cbind(Xlo[-1], Xhi[-1])
    # JNlohigrps <- cbind(Xlo[-1], Xhi[-1])
    # JNlohigrps <- cbind(XL1, XL2)
    colnames(JNlohigrps) <- c('Low Value', 'High Value')
    rownames(JNlohigrps) <- compnoms[-1,]
    
    if (verbose) {
      message('\n\nSimultaneous regions of significance -- Johnson-Neyman Technique:')
      message('\nUsing Bonferroni F, p = .05, 2-tailed; see Huitema, 1980, p. 293')
      message('The regression lines for group comparisons are significantly different')
      message('at IDV scores < Low Value & > High Value:\n')
      
      print(round(JNlohigrps,3), print.gap = 4)
      if (any(is.na(JNlohigrps))) message('NA indicates that meaningful values could not be computed.\n')	
    }
  }		
  
  
  
  
  
  ######################################    plot data    #######################################
  
  
  plotdon <- NULL
  
  if (plot_type == 'interaction') {
    
    
    # IV min & max values for plot -- categorical MOD
    if (MOD_type == 'factor') {	
      
      Ngroups <- length(levels(donnes[,MOD]))
      plotdon <- rep(-9999,Ngroups)
      
      # IV min & max values for plot -- dichotomous IV
      if (IV_type == 'factor' & length(levels(donnes[,IV])) == 2) {
        IV_min <- rep( (levels(donnes[,MOD])[1] = 0), Ngroups)
        IV_max <- rep( (levels(donnes[,MOD])[2] = 1), Ngroups)
      }
      
      # IV min & max values for plot -- numeric IV      
      if (IV_type == 'numeric') {
        IV_min <- IV_max <- NULL
        for (lupeD in 1:Ngroups) {
          
          dontemp <- subset(donnes, mod_values==levels(donnes[,MOD])[lupeD], select=IV)
          
          if (IV_range == 'quantiles' | IV_range == 'tumble') {   # using the 10th & 90th
            IVquants <- quantile(dontemp, na.rm=T, probs=quantiles_IV)	
            IV_min <- c(IV_min, IVquants[1])
            IV_max <- c(IV_max, IVquants[2])
          } else if (IV_range == 'minmax') {  
            IV_min <- c(IV_min, min(dontemp))
            IV_max <- c(IV_max, max(dontemp))
          } else if (IV_range == 'AikenWest') { 
            IVmn <- sapply(dontemp[IV], mean, na.rm = TRUE)
            IVsd <- sapply(dontemp[IV], sd, na.rm = TRUE)
            IV_min <- c(IV_min, (IVmn - IVsd))
            IV_max <- c(IV_max, (IVmn + IVsd))
          } else if (IV_range == 'numeric') {	
            IV_min <- c(IV_min, IV_range_user[1])
            IV_max <- c(IV_max, IV_range_user[2])
          }			
        }		
      }
      
      # the dummy codes for MOD
      MODdummys <- contrasts(donnes[,MOD])
      
      plotdon <- rbind( cbind(IV_min, MODdummys), cbind(IV_max, MODdummys) )
      colnames(plotdon) <- c(IV, MODnew)
    }
    
    
    
    # IV min & max values for plot -- continuous MOD
    if (MOD_type == 'numeric') {	
      
      # IV min & max values for plot -- dichotomous IV
      if (IV_type == 'factor' & length(levels(donnes[,IV])) == 2) {
        IV_min <- levels(donnes[,MOD])[1] = 0
        IV_max <- levels(donnes[,MOD])[2] = 1
      }
      
      # IV min & max values for plot -- numeric IV
      if (IV_type == 'numeric') {
        
        plotdon <- rep(-9999,2)
        if (IV_range == 'tumble') { # | 
          # Bodner 2016 p 598 tumble graph method for IV ranges
          # Use Equation 3 to predict the conditional mean values of the target predictor 
          # X for each of the moderator variable values chosen in Step 1. The square root 
          # of the mean square residual from the analysis of variance summary table for 
          # this model is an estimate of the SD of the residuals around the predicted values; 
          # this value is added and subtracted from each predicted target variable value 
          # resulting in two values of the target variable X for each chosen moderator variable value. 
          formB <- as.formula(paste(IV, MOD, sep=" ~ "))	
          eq3 <- lm(formB, donnes)
          sumtab <- summary(eq3)
          sqrmsr <- sumtab$sigma	
          for (lupe in 1:length(mod_values)) {
            predval <- sumtab$coeff[1,1] + sumtab$coeff[2,1] * mod_values[lupe]
            IV_min <- predval - sqrmsr
            plotdon <- rbind( plotdon, c(IV_min, mod_values[lupe]))
            IV_max <- predval + sqrmsr
            plotdon <- rbind( plotdon, c(IV_max, mod_values[lupe]))		
          }
          plotdon <- plotdon[-1,]		
        } else if (IV_range == 'quantiles') { # using the 10th & 90th
          IVquants <- quantile(donnes[IV], na.rm=T, probs=quantiles_IV)	
          IV_min <- IVquants[1]
          IV_max <- IVquants[2]
        } else if (IV_range == 'minmax') { 
          IV_min <- min(donnes[IV])
          IV_max <- max(donnes[IV])
        } else if (IV_range == 'AikenWest') { 
          IVmn <- sapply(donnes[IV], mean, na.rm = TRUE)
          IVsd <- sapply(donnes[IV], sd, na.rm = TRUE)
          IV_min <- IVmn - IVsd
          IV_max <- IVmn + IVsd
        } else if (IV_range == 'numeric') {	
          IV_min <- IV_range_user[1]
          IV_max <- IV_range_user[2]
        }
        
        if (min(plotdon) == -9999) plotdon <- expand.grid(c(IV_min,IV_max), mod_values)
        
        colnames(plotdon) <- c(IV, MODnew)
      }
    }
    
    
    
    # add COVARS data to plotdon for the predict function -- using the means of each
    if (!is.null(COVARS)) {
      # COVARSmn <- sapply(donnes[COVARS], mean, na.rm = TRUE)
      
      COVARSmn <- matrix(sapply(donnes[COVARS], mean, na.rm = TRUE), nrow=1)
      
      COVARSmn <- matrix(rep(COVARSmn,nrow(plotdon)), nrow=nrow(plotdon), ncol=4, byrow=TRUE)
      plotdon <- cbind(plotdon, COVARSmn)
      colnames(plotdon)[3:(2+length(COVARS))] <- COVARS
      
      # plotdon <- cbind(plotdon, rep(COVARSmn,nrow(plotdon)))
      # colnames(plotdon)[3:(2+length(COVARS))] <- COVARS
    }
    
    plotdon <- data.frame(plotdon)
    predvals <- predict(modelXN, type='response', plotdon)
    # removing COVARS from plotdon because not needed for plot
    plotdon <- subset(plotdon, select=c(IV,MODnew))
    # adding the predicted DV values
    plotdon$predDV <- predvals
    colnames(plotdon)[colnames(plotdon) == 'predDV'] <- DV
    
    if (MOD_type == 'factor') { 
      plotdon <- plotdon[,c(IV,DV)]
      plotdon[,MOD] <- rep(mod_values,2) 
    }
    
    
    # set the range for the x and y axis 
    xrange <- range(plotdon[IV]) 
    
    yrange <- range(plotdon[DV])
    if (!is.null(DV_range))  yrange <- DV_range 
    
    
    # set up the plot 
    
    if (is.null(Xaxis_label))   Xaxis_label <- IV
    
    if (is.null(Yaxis_label))   Yaxis_label <- DV
    
    if (is.null(plot_title))    plot_title <- 'Interaction Plot'
    
    if (is.null(legend_label))  legend_label <- MOD
    
    
    plot(xrange, yrange, type="n", xlab=Xaxis_label, ylab=Yaxis_label, cex.lab=1.3, main=plot_title ) 
    
    for (i in 1:length(mod_values)) {
      dum <- subset(plotdon, plotdon[,MOD]==mod_values[i], select = c(IV,DV))
      lines(dum, type="b", lwd=1.5, lty=1, col=i, pch=19); #points(dum)
    }
    
    if (MOD_type == 'numeric') {legvalues <- round(mod_values,2)} else {legvalues <- mod_values}
    legend("topleft", legend=legvalues, title=legend_label, col=1:length(mod_values), 
           bty="n", lty=1, lwd=2) # ,inset = c(.60,.03)
    
  } # end of  if (plot_type == 'interaction') 
  
  
  
  output <- list(modelMAINsum=modelMAINsum, mainRcoefs=mainRcoefs, 
                 modeldata=modeldata, collin_diags=collin_diags,
                 modelXN=modelXN, modelXNsum=modelXNsum, 
                 RsqchXn=RsqchXn, fsquaredXN=fsquaredXN, xnRcoefs=xnRcoefs, 
                 simslop=simslop, simslopZ=simslopZ, 
                 plotdon=plotdon, JN.data = JN.data, ros=ros, DV=DV, IV=IV, MOD=MOD, family='OLS')
  class(output) <- "MODERATED_REGRESSION"
  
  
  
  # plot of Johnson-Neyman regions of significance
  if (plot_type == 'regions') {
    
    if (is.null(Xaxis_label))   Xaxis_label <- MOD
    
    if (is.null(Yaxis_label))   Yaxis_label <- c(paste("Simple Slopes of",IV,'on',DV))
    
    if (is.null(plot_title))    plot_title <- c("Simple Slopes of",IV,'on',paste(DV,' by ',MOD,sep=""))
    
    if (is.null(legend_label))  legend_label <- 'Simple Slope'
    
    REGIONS_OF_SIGNIFICANCE(model=output,  
                            IV_range=IV_range, MOD_range=MOD_range, 
                            plot_title=plot_title, Yaxis_label=Yaxis_label, 
                            Xaxis_label=Xaxis_label, legend_label=legend_label,
                            names_IV_MOD=NULL) 
  }
  
  if (plot_type == 'residuals') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(16,2,3,4))
  
  if (plot_type == 'diagnostics') 
    diagnostics_plots(modelMAIN=modelMAIN, modeldata=modeldata, plot_diags_nums=c(9,12,13,14))
  
  
  return(invisible(output))
  
}



