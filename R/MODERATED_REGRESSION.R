

MODERATED_REGRESSION <- function (data, DV, IV, MOD,
                                  IV_type = 'numeric', IV_range='tumble',
                                  MOD_type = 'numeric', MOD_levels='quantiles', 
                                  MOD_range=NULL, MOD_reflevel=NULL,
                                  quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                                  COVARS = NULL,
                                  center = TRUE,  
                                  CI_level = 95,
                                  MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                                      thin = 1, burnin = 1000, 
                                                      HDI_plot_est_type = 'raw'),
                                  plot_type = 'residuals', plot_title=NULL, DV_range=NULL,
                                  Xaxis_label=NULL, Yaxis_label=NULL, legend_label=NULL,
                                  JN_type = 'Huitema',
                                  verbose=TRUE ) {
  
  
  "<-<-" <- NULL   # need this or else get "no visible global function definition for '<-<-' " on R CMD check
  
  if ( is.null(IV) | is.null(MOD) ) {
    message('\n\nThe data for the analyses were not properly specified.\n')
    message('\nThe IV or MOD arguments were NULL. Expect errors.\n')
  }
  
  donnes <- data[,c(DV,IV,MOD,COVARS)]
  
  # if IV_type = 'factor' but the IV data are not a factor, then convert it to a factor
  if (IV_type == 'factor' & !is.factor(donnes[,IV]))  donnes[,IV] <- factor(donnes[,IV])
  
  # if MOD_type = 'factor' but the MOD data are not a factor, then convert it to a factor
  if (MOD_type == 'factor' & !is.factor(donnes[,MOD]))  donnes[,MOD] <- factor(donnes[,MOD])
  
  
  # if IV_type = 'numeric' or 'integer' but the IV data are not numeric or integer, then convert it to a factor
  if ( (IV_type == 'numeric' | IV_type == 'integer') &
       (!is.numeric(donnes[,IV]) & !is.integer(donnes[,IV])) )  {
    donnes[,IV] <- factor(donnes[,IV])
    IV_type <- 'factor'
    message('\nIV_type was specified as numeric but the IV data are not numeric, so it was converted to a factor')
  }
  
  # if MOD_type = 'numeric' or 'integer' but the MOD data are not numeric or integer, then convert it to a factor
  if ( (MOD_type == 'numeric' | MOD_type == 'integer') &
       (!is.numeric(donnes[,MOD]) & !is.integer(donnes[,MOD])) )  {
    donnes[,MOD] <- factor(donnes[,MOD])
    MOD_type <- 'factor'
    message('\nMOD_type was specified as numeric but the MOD data are not numeric, so it was converted to a factor')
  }
  
  
  # clean up IV, MOD, COVARS factors, contrasts
  IV_type_cleanup(donnes, allIVnoms = c(IV, MOD, COVARS)) 
  
  
  if (is.numeric(IV_range)) {
    IV_range_user <- IV_range
    IV_range <- 'numeric'
  }
  
  if (anyNA(donnes)) {donnes <- na.omit(donnes); NAflag = TRUE} else {NAflag = FALSE}	
  
  if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
  
  
  if (length(MOD_levels) > 1) {
    
    mod_values <- MOD_levels
    
    MOD_levels <- 'user'
  }
  
  
  # centering, if requested
  if (center) {
    if (is.numeric(donnes[,IV]))   donnes[,IV]  <- donnes[,IV]  - mean(donnes[,IV])	
    if (is.numeric(donnes[,MOD]))  donnes[,MOD] <- donnes[,MOD] - mean(donnes[,MOD])	
  }
  
  
  # if MOD is a factor, set the baseline
  if (is.factor(donnes[,MOD])) {
    
    # For factors in R, the baseline group is normally the first-appearing value  
    # in the df. Instead, in the following code,
    # (1) the MOD_reflevel, if provided, will be used as the baseline; or
    # (2) if the terms "control" or "Control" or "baseline" or "Baseline" 
    #     appear in the names of a factor level, then that factor level will be
    #     used as the baseline; 
    # (3) otherwise the baseline will be the earliest of the alphabetically-
    #     ordered factor levels
    
    fac_levels <- levels(donnes[,MOD])
    
    if (!is.null(MOD_reflevel))  {
      # test if MOD_reflevel is in fact one of the fac_levels
      if ((MOD_reflevel %in% fac_levels) == FALSE) {
        message('\n\nThe name provided on the MOD_reflevel argument is not one of the')
        message('MOD factor levels. The default reference level will be used instead.\n')
        MOD_reflevel <- NULL
      } else {base_level <- MOD_reflevel}
    }
    
    if (is.null(MOD_reflevel)) {
      if (length(grep(pattern = 'Control|control|baseline|Baseline', x = fac_levels)) > 0) {
        base_level <- grep(pattern = 'Control|control|baseline|Baseline', x = fac_levels)
      } else {
        base_level <- sort(fac_levels)[1]
      }
    }
    
    # levels(donnes[,MOD])
    donnes[,MOD] <- relevel(donnes[,MOD], ref = base_level)
    # levels(donnes[,MOD])
  }
  
  # descriptives
  allIVnoms <- c(DV, IV, MOD)
  if (!is.null(COVARS))  allIVnoms <- c(allIVnoms, COVARS)
  if (verbose)  descriptives(donnes, allIVnoms)
  
  
  for (lupe in 1:2) {
    
    if (lupe == 1)  prevRsq <- prevModel <- prevmod_lmBF <- NULL
    
    # keeping info on previous model, if there was one, for Rsq change stats
    if (lupe > 1) {
      prevModel <- xx$model
      prevRsq   <- xx$Rsq
      if (MCMC_options$MCMC) prevmod_lmBF <- xx$mod_lmBF
    }
    
    message('\n\nBlock ', lupe, ':')	
    
    if (lupe == 1) {
      
      formule <- as.formula(paste(DV, paste(c(IV,MOD,COVARS), collapse=" + "), sep=" ~ "))
      
      # noms_list <- list(DV=DV, IV=IV, MOD=MOD, COVARS=COVARS)
      
      xx <- reg_stats_boc(formule = formule, data = donnes, CI_level = CI_level, 
                          prevRsq = NULL, prevModel = NULL, prevmod_lmBF = NULL,
                          MCMC_options = MCMC_options) #, noms_list=NULL) #noms_list)

      if (verbose) {
        
        message('\n\nThe DV is: ', DV)
        message('\nThe IV is: ', IV)
        message('\nThe MOD is: ', MOD)
        if (!is.null(COVARS))  message('\nThe covariates are: ', COVARS)
        message('\n\nMultiple regression statistics for when the product term is not in the model:')
        
        resultats(xx)
      }
    }
    
    if (lupe == 2) {
      
      prevModel <- xx$model
      prevRsq   <- xx$Rsq
      if (MCMC_options$MCMC) prevmod_lmBF <- xx$mod_lmBF

      # the term names for IV 
      if (IV_type != 'factor')  IV_terms <- IV
      if (IV_type == 'factor')  {
        formC <- as.formula(paste('~', IV, sep=" "))
        IV_terms <- colnames(model.matrix(formC, data = donnes))[-1]
      }
      
      # the term names for MOD 
      if (MOD_type != 'factor')  MOD_terms <- MOD
      if (MOD_type == 'factor')  {
        formC <- as.formula(paste('~', MOD, sep=" "))
        MOD_terms <- colnames(model.matrix(formC, data = donnes))[-1]
      }
      
      # the term names for the COVARS 
      COVARS_terms <- NULL
      if (length(COVARS) > 0) {
        for (lupe in 1:length(COVARS)) {
          if (is.numeric(donnes[,COVARS[lupe]]) | 
              is.integer(donnes[,COVARS[lupe]])) {
            COVARS_terms <- c(COVARS_terms, COVARS[lupe])
          } else { 
            formC <- as.formula(paste('~', COVARS[lupe], sep=" "))
            COVARS_terms <- c(COVARS_terms, colnames(model.matrix(formC, data = donnes))[-1]) }
        }
      }
      
      if (is.null(COVARS))
        formXn <- as.formula(paste(DV, paste(IV, ' * ', MOD, sep = ""), sep=" ~ "))
      
      if (!is.null(COVARS))
        formXn <- as.formula(paste(DV,   paste(
          paste(IV, ' * ', MOD, sep = ""),
          paste(COVARS, collapse=" + ", sep = ""), 
          sep=' +'), sep=" ~ "))
      
      xx <- reg_stats_boc(formule = formXn, data = donnes, CI_level = CI_level, 
                          prevRsq = NULL, prevModel = NULL, prevmod_lmBF = NULL,
                          MCMC_options = MCMC_options) #, noms_list=NULL) #noms_list)
      
      PROD_terms <- setdiff( rownames(xx$partialRcoefs), c(IV, MOD_terms, COVARS_terms))
      
      # noms_list <- list(DV=DV, IV=IV, MOD_terms=MOD_terms, PROD_terms=PROD_terms, COVARS=COVARS)
      
      if (verbose) {
        
        message('\n\nThe DV is: ', DV)
        message('\nThe IV is: ', IV)
        message('\nThe MOD is: ', MOD)
        if (!is.null(COVARS))  message('\nThe covariates are: ', COVARS)
        message('\n\nMultiple regression statistics for when the product term(s) is in the model:')
        resultats(xx)
      }  
    }
  }
  
 
  
  
  # the values/levels of MOD, if it is continuous
  if (MOD_type != 'factor') {  
    
    if (MOD_levels=='quantiles')  
      mod_values <- quantile(donnes[,MOD], probs=quantiles_MOD)
    
    if (MOD_levels=='AikenWest') {
      MODmn <- sapply(donnes[MOD], mean, na.rm = TRUE)
      MODsd <- sapply(donnes[MOD], sd, na.rm = TRUE)
      MODlo <- MODmn - MODsd
      MODhi <- MODmn + MODsd
      mod_values <- c(MODlo, MODmn, MODhi)
    }	
    
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
  
  
  
  modelXN <- xx$model
  modeldata <- xx$modeldata
  

  ###############################  simple slopes  ##################################
  
  simslop <- simslopZ <- NULL
  
  # only if IV is not a factor
  if (IV_type == 'numeric') {
    
    coefs <- modelXN$coefficients
    Sb <- vcov(modelXN)[2:length(coefs),2:length(coefs)]
    if (MOD_type != 'factor') {
      slopes   <- coefs[IV] + coefs[PROD_terms] * mod_values
      intercepts  <- coefs[MOD] * mod_values + coefs['(Intercept)']
      SEslopes <- sqrt( Sb[IV,IV] + 2*mod_values*Sb[IV,PROD_terms] +  mod_values**2 * Sb[PROD_terms,PROD_terms])
    } else { 
      slopes   <- c(coefs[IV], (coefs[IV] + coefs[PROD_terms]))
      intercepts  <- c(coefs[1],  (coefs[1]  + coefs[MOD_terms]))
      diagSB   <- diag(Sb)
      SEslopes <- c(sqrt(Sb[1,1]), sqrt( abs(Sb[1,1] - diagSB[PROD_terms])))
    }
    tslopes <- slopes / diag(SEslopes)
    tslopes <- slopes / SEslopes
    N <- nrow(donnes)
    df <- modelXN$df.residual
    k <- length(coefs) - 1     # number of predictors 
    dfs <-  N - k -1  
    pslopes <- (1 - pt(abs(tslopes),dfs)) * 2
    # CIs - 2003 Cohen Aiken West p 278
    tabledT <- qt(((100 - CI_level) / 2 * .01), dfs)   # tabledT <- qt(.975, dfs)  
    me <- tabledT * SEslopes
    confidLo <- slopes - me
    confidHi <- slopes + me
    
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
        split(data.frame(donnes[,DV], donnes[,IV]), donnes[,MOD]),  
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
  }
  
  
  
  ###########################  Johnson-Neyman regions of significance  ###########################
  
  JN.data <- ros <- NULL
  
  # only if IV is not a factor
  if (IV_type == 'numeric') {
    
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
    
    
    B <- modelXN$coefficients[c('(Intercept)',IV,MOD,PROD_terms)]
    S <- stats::vcov(modelXN)[c('(Intercept)',IV,MOD,PROD_terms), c('(Intercept)',IV,MOD,PROD_terms)]
    
    tcrit <- qt(((100 - CI_level) / 2 * .01), modelXN$df.residual)  # tcrit <- qt(0.975, modelXN$df.residual)
    
    # construct the quadratic equation
    a <- tcrit^2 * S[PROD_terms,PROD_terms] - B[PROD_terms]^2
    b <- 2 * (tcrit^2 * S[IV,PROD_terms] - B[IV]*B[PROD_terms])
    c <- tcrit^2 * S[IV,IV] - B[IV]^2
    
    
    # when JN values can be computed
    if (b^2-4*a*c >= 0) {
      JNa <- (-b - sqrt(b^2-4*a*c)) / (2*a)
      JNb <- (-b + sqrt(b^2-4*a*c)) / (2*a)
      ros <- sort(c(JNa,JNb))
      
      # adding the ros values to MODvalues, but only if they are within the MOD min/max range
      MODvalues <- sort(c(MODvalues, ros[ros >= MOD_min & ros <= MOD_max ])	)
      
      # data for plot
      SimpleSlope <- B[IV]+(B[PROD_terms]*MODvalues)
      StdError <- sqrt((S[IV,IV])+(MODvalues^2*S[PROD_terms,PROD_terms])+(2*MODvalues*S[IV,PROD_terms]) )   
      CI.L <- SimpleSlope - tcrit*StdError  
      CI.U <- SimpleSlope + tcrit*StdError  
      JN.data <- data.frame(SimpleSlope, CI.L, CI.U, MODvalues, StdError)
      
      if (verbose) {
        
        message('\n\nJohnson-Neyman regions of significance interval:')  
        message('\n   The slopes for ',IV,' predicting ',DV,' are significant when ',MOD)
        message('\n','   values are outside of (lower and higher than) this range: ',round(ros[1],2),' to ',round(ros[2],2))
        
        # slope & CI values at the ros points
        if (ros[1] >= MOD_min) {
          # if ros[1] is > the largest JN.data$MODvalues, then use the largest
          if (ros[1] > tail(JN.data$MODvalues, 1)) { rownum <- length(JN.data$MODvalues)
          } else {rownum <- which(JN.data$MODvalues == ros[1])} # the row with the high ros MOD value info
          message('\n   When ',MOD,' is ',round(JN.data$MODvalues[rownum],2),
                  ':   slope = ',round(JN.data$SimpleSlope[rownum],2),
                  '  CI_lb = ',round(JN.data$CI.L[rownum],2),
                  '  CI_ub = ',round(JN.data$CI.U[rownum],2),
                  '  SE = ',round(JN.data$StdError[rownum],2))
        }
        if (ros[2] <= MOD_max) {
          # if ros[2] is < the smallest JN.data$MODvalues, then use the smallest
          if (ros[2] < JN.data$MODvalues[1]) { rownum <- 1
          } else {rownum <- which(JN.data$MODvalues == ros[2])} # the row with the high ros MOD value info
          message('\n   When ',MOD,' is ',round(JN.data$MODvalues[rownum],2),
                  ':    slope = ',round(JN.data$SimpleSlope[rownum],2),
                  '    CI_lb = ',round(JN.data$CI.L[rownum],2),
                  '   CI_ub = ',round(JN.data$CI.U[rownum],2),
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
    Ngroups <- length(MOD_terms) + 1  # OK?	 	
    Ncombins <- choose(Ngroups,2)  # the # of possible 2-group comparisons
    grpNs <- grp_IV_MNs <- grp_DV_MNs <- sscp11 <- sscp12 <- sscp22 <- ssreg <- 
      intercepts <- slopes <- ssresid <- rep(-9999, Ngroups)
    JNlohigrps <- c(NA,NA)
    compnoms <- NA
    # Xhi <- Xlo <- NA  # matrix(NA, choose(Ngroups,2), 2) 
    alpha <- ((100 - CI_level) / 2 * .01)  # .05
    
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
        
        if ((B**2 - A*C) > 0) {
          
          hi <- (-B + (sqrt(B**2 - A*C)) ) / A
          lo <- (-B - (sqrt(B**2 - A*C)) ) / A
          
          if (hi > lo) { JNlohigrps <- rbind(JNlohigrps, c(lo,hi))
          } else { JNlohigrps <- rbind(JNlohigrps, c(NA,NA)) }
        } else { JNlohigrps <- rbind(JNlohigrps, c(NA,NA)) }
        
        compnoms <- rbind(compnoms,  paste('Groups',lupe1,'&',lupe2, sep=' ') )
      }
    }
    
    JNlohigrps <- JNlohigrps[-1, ,drop = FALSE]
    colnames(JNlohigrps) <- c('Low Value', 'High Value')
    rownames(JNlohigrps) <- compnoms[-1,]
    

    
    if (verbose) {
      message('\n\nSimultaneous regions of significance -- Johnson-Neyman Technique:')
      message('\nUsing Bonferroni F, p = .05, 2-tailed; see Huitema, 1980, p. 293')
      message('The regression lines for group comparisons are significantly different')
      message('at IDV scores < Low Value & > High Value:\n')
      
      print(round(JNlohigrps,3), print.gap = 4)
      if (any(is.na(JNlohigrps))) message('\n    NA indicates that meaningful values could not be computed.\n')	
    }
  }		
}
  
  
  ###########################   interaction plot    ############################
  
  plotdon <- NULL
  
  if (plot_type == 'interaction') {
    
    # plotdon will have more IV values when IV_range == 'tumble' than for the other
    # IV_range methods. Therefore, plotdon will be computed in different ways:
    
    # can't compute tumble graph ranges when IV is a factor
    if (IV_range == 'tumble' & IV_type == 'factor')  IV_range = 'AikenWest'
    
    if (IV_range == 'tumble') {
      
      if (MOD_type == 'factor') {
        
        # Bodner 2016 p 598 tumble graph method for IV ranges
        # For categorical moderator variables, researchers will likely have a distribu-
        # tion of scores on the target predictor X within each category of M, and this
        # feature provides several options. For example, the two values of X can be
        # selected using the same quantities computed separately within each category
        # (e.g., using the first and third quartile values in a category or at +1 within-
        # category SD from the category mean for X). If the target variable variance is
        # reasonably homogeneous across levels of the categorical moderator variable
        # (e.g., the ratio of the largest to the smallest variance being less than 2), the SD
        # can be based on the square root of the pooled (within-moderator level) target
        # variable variances.
        
        # I am using the withing group mean & +1 or -1 SD method, with the SD being
        # the pooled SD (from lm)
        
        formB <- as.formula(paste(IV, MOD, sep=" ~ "))
        eq3 <- lm(formB, donnes)
        sumtab <- summary(eq3)
        sqrmsr <- sumtab$sigma  # the pooled SD
        IV_vals <- MOD_vals <- NULL
        
        for (lupe in 1:length(mod_values)) {   # for (lupeD in 1:Ngroups) {
          dontemp <- subset(donnes, mod_values==levels(donnes[,MOD])[lupe], select=IV)
          IVmn <- sapply(dontemp[IV], mean, na.rm = TRUE)
          IV_min <- IVmn - sqrmsr
          IV_max <- IVmn + sqrmsr
          IV_vals <- c(IV_vals, c(IV_min, IV_max))
          MOD_vals <- c(MOD_vals, rep(mod_values[lupe],2))
        }
        plotdon <- data.frame(IV = IV_vals, MOD = MOD_vals)
        colnames(plotdon) <- c(IV, MOD)
      }
      
      if (MOD_type == 'numeric') {
        
        # Bodner 2016 p 598 tumble graph method for IV ranges
        # For continuous moderator variables, the options are more limited, given that
        # there may be no or few data values of X at a given value of M. If the relationship
        # between X and M is approximately linear and the residuals of X about the
        # regression line of X on M are reasonably symmetric and homoscedastic, then
        # researchers may estimate quantities in the conditional distribution of X given M
        # (e.g., at +1 residual SD around the predicted value of X given M). In particular,
        # regress the target predictor X on the moderator variable M to obtain a regression
        # equation, that is,
        # Equation 3
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
        plotdon <- NULL
        IV_vals <- MOD_vals <- NULL
        
        for (lupe in 1:length(mod_values)) {
          predval <- sumtab$coeff[1,1] + sumtab$coeff[2,1] * mod_values[lupe]
          IV_min <- predval - sqrmsr
          IV_max <- predval + sqrmsr
          IV_vals <- c(IV_vals, c(IV_min, IV_max))
          MOD_vals <- c(MOD_vals, rep(mod_values[lupe],2))
        }
        plotdon <- data.frame(IV = IV_vals, MOD = MOD_vals)
        colnames(plotdon) <- c(IV, MOD)
      }
    }
    
    if (IV_range != 'tumble') {
      
      # IV values for plotdon
      if (IV_type == 'factor')  IVvalues <- levels(donnes[,IV])
      
      if (IV_type == 'numeric') {
        
        if (IV_range == 'quantiles') { 
          IVquants <- quantile(donnes[IV], na.rm=T, probs=quantiles_IV)	
          IV_min <- IVquants[1]
          IV_max <- IVquants[2]
        }
        
        if (IV_range == 'minmax') { 
          IV_min <- min(donnes[IV])
          IV_max <- max(donnes[IV])
        }
        
        if (IV_range == 'AikenWest') { 
          IVmn <- sapply(donnes[IV], mean, na.rm = TRUE)
          IVsd <- sapply(donnes[IV], sd, na.rm = TRUE)
          IV_min <- IVmn - IVsd
          IV_max <- IVmn + IVsd
        }
        
        if (IV_range == 'numeric') {	
          IV_min <- IV_range_user[1]
          IV_max <- IV_range_user[2]
        }
        
        IVvalues <- c(IV_min, IV_max)
      }  # end of  if (IV_type == 'numeric')
      
      
      # MOD values for plotdon = "mod_values", computed above
      
      plotdon <- expand.grid(IVvalues, mod_values)
      colnames(plotdon) <- c(IV, MOD)
    }

    
    # add COVARS data to plotdon for the predict function
    if (!is.null(COVARS)) {
      
      # COVARSmn <- matrix(sapply(donnes[COVARS], mean, na.rm = TRUE), nrow=1)
      # 
      # COVARSmn <- matrix(rep(COVARSmn,nrow(plotdon)), nrow=nrow(plotdon), 
      #                    ncol=length(COVARS), byrow=TRUE)
      # colnames(COVARSmn) <- COVARS
      # 
      # plotdon <- cbind(plotdon, COVARSmn)
      # # colnames(plotdon)[3:(2+length(COVARS))] <- COVARS
      
      # COVARS_values_list <- c()
      
      # when there are factors, cycle through the COVARS, use mean if continuous, use baseline if categorical
      if (length(COVARS) > 0) {
        
        factor_vars <- names(xx$model$xlevels)
        
        for (lupe in 1:length(COVARS)) {
          
          if (COVARS[lupe] %in% factor_vars) {
            
            # COVARS_values_list[[COVARS[lupe]]] <- sort(list_xlevels[[COVARS[lupe]]])[1]
            
            # COVARS_values_list[[COVARS[lupe]]] <- unlist(modelXN$xlevels)[1]
            
            plotdon[,COVARS[lupe]] <- unlist(modelXN$xlevels)[1]
            
          } else { 
            
            # COVARS_values_list[[COVARS[lupe]]] <- mean(xx$modeldata[,COVARS[lupe]]) 
            
            # plotdon$grade90_factor <- unlist(modelXN$xlevels)[1] }
            
            plotdon[,COVARS[lupe]] <- mean(xx$modeldata[,COVARS[lupe]])
          }
        }
      }
    }
    
    
    predvals <- predict(modelXN, type='response', plotdon)
    # removing COVARS from plotdon because not needed for plot
    # plotdon <- subset(plotdon, select=c(IV,MOD_terms))
    plotdon <- subset(plotdon, select=c(IV,MOD))
    # adding the predicted DV values
    plotdon$predDV <- predvals
    colnames(plotdon)[colnames(plotdon) == 'predDV'] <- DV
    
    # set the range for the x and y axis 
    if (IV_type == 'numeric') xrange <- range(plotdon[IV]) 
    
    yrange <- round(range(plotdon[DV]), 2)
    # increase the upper ylim by 50%
    if (IV_type == 'factor')  yrange <- better_ylim(yrange, buffer = 0.45)
    
    if (!is.null(DV_range))  yrange <- DV_range 
    
    # set up the plot 
    
    if (is.null(Xaxis_label))   Xaxis_label <- IV
    
    if (is.null(Yaxis_label))   Yaxis_label <- DV
    
    if (is.null(plot_title))    plot_title <- 'Interaction Plot'
    
    if (is.null(legend_label))  legend_label <- MOD
    
    
    plotfun(testdata=plotdon, list_xlevels = modelXN$xlevels, 
            DV_predicted=DV, CIs=FALSE, kind='OLS',
            xlim=xrange, ylim=yrange, xlab=Xaxis_label, ylab=Yaxis_label, 
            title=plot_title, cols_user=NULL,
            IV_focal_1=IV, IV_focal_1_values=IVvalues, #plotdon[,IV], 
            IV_focal_2=MOD, IV_focal_2_values=mod_values) 

  } # end of  if (plot_type == 'interaction') 
  
  
  
  output <- list(modelsum=xx$modelsum, partialRcoefs=xx$partialRcoefs, 
                 modeldata=xx$modeldata, collin_diags=xx$collin_diags,
                 model=xx$model, modelsum=xx$modelsum, 
                 Rsqch=xx$Rsqch, fsquared=xx$fsquared, partialRcoefs=xx$partialRcoefs, 
                 simslop=simslop, simslopZ=simslopZ, 
                 plotdon=plotdon, JN.data = JN.data, ros=ros, DV=DV, IV=IV, MOD=MOD, 
                 family='OLS', # noms_list=noms_list, 
                 chain_dat = xx$chains, 
                 Bayes_HDIs = xx$Bayes_HDIs)
  
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
    diagnostics_plots(model=xx$model, modeldata=modeldata, plot_diags_nums=c(16,2,3,4))
  
  if (plot_type == 'diagnostics') 
    diagnostics_plots(model=xx$model, modeldata=modeldata, plot_diags_nums=c(9,12,13,14))
  
  if (plot_type == 'Bayes_HDI' & MCMC_options$MCMC & !is.null(xx$chains)) 
    Bayes_HDI_plot(chain_dat = xx$chains, # [,c(IV, MOD_terms, PROD_terms)] 
                   Bayes_HDIs = xx$Bayes_HDIs,
                   CI_level = CI_level, # SDs = xx$SDs[c(DV, IV, MOD_terms, PROD_terms)],
                   HDI_plot_est_type = MCMC_options$HDI_plot_est_type)
    
  return(invisible(output))
  
}



