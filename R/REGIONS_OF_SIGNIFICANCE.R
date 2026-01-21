


# J-N for lm & lme (nlme package) interaction models

REGIONS_OF_SIGNIFICANCE <- function(model,  
                                    IV_range=NULL, MOD_range=NULL, 
                                    plot_title=NULL, Xaxis_label=NULL, Yaxis_label=NULL, 
                                    legend_label=NULL,
                                    names_IV_MOD=NULL) { 
  
  
  # COVARS ??? FIX
  
  # JN_type = 'Huitema',
  
  
  # MOD must be numeric
  
  # If IV & MOD are the only predictors, in that order, then there is no need to
  # specify names_IV_MOD_raw. The J-N computations will be based on these variables being
  # in these positions.
  
  # If there are covariates in addition to IV & MOD as predictors, or if the order
  # is not IV then MOD, then names_IV_MOD_raw must be specified.
  
  # The names_IV_MOD_model argument can be used to id the key terms () from an lme model
  # that involved more than IV, MOD, & Xn terms. The argument is used only to create
  # the key B and S objects for the J-N analyses. Other terms in the model are ignored.
  
  # If model is an lme object & IV is a two-level factor, then names_IV_MOD_model
  # must be specified (because lme alters the variable names).
  
  
  
  if (inherits(model, 'MODERATED_REGRESSION')) {    # if (class(model) == 'MODERATED_REGRESSION') {
    
    DV  <- model$DV
    IV  <- model$IV
    MOD <- model$MOD
    
    JN.data <- model$JN.data	
    
    ros <- model$ros	
    
    if (is.null(JN.data)) {
      message('\n\nA MODERATED_REGRESSION model was provided, but it does not contain data')
      message('\nfor a Johnson-Neyman regions of significance plot.')
      message('\nA possible reason is that the MOD (moderator) variable was not continuous.')	
      message('\nA Johnson-Neyman regions of significance plot is provided only for')
      message('\ncontinuous/numeric moderators.')	
    }
    
  }
  
  
  
  if (inherits(model, 'lme')) {   # if (class(model) == 'lme') {
    
    rawdata <- getData(model)
    
    noms <- names(fixef(model))
    
    DV <- rownames(attributes(terms(model))$factors)[1]
    
    # if names_IV_MOD = NULL, the IV & MOD names are taken from the lme model, from the 2nd & 3rd terms
    if (is.null(names_IV_MOD)) { 
      IV  <- noms[2]
      MOD <- noms[3]
      
      # IV  <- rownames(attributes(terms(model))$factors)[2]
      # MOD <- rownames(attributes(terms(model))$factors)[3]
    } else {
      IV  <- names_IV_MOD[1]
      MOD <- names_IV_MOD[2]		
    }
    
    # variable/column names issue, now bypassed: The names from names(fixef(model)) are the ones
    # from lme, & lme sometimes changes the names, e.g., for factors -- they are 
    # not always the same as the names in rawdata
    
    # find the name of the product term - that appears in the lme model 
    # lme changes the orig data names of factors, in which case the model names != the data names
    noms <- names(fixef(model))
    PRODpos <- which(grepl(IV, noms) & grepl(MOD, noms) & grepl(':', noms) )
    PROD <- noms[PRODpos]
    
    # IV min & max 
    if (is.null(IV_range)) {
      
      # determine if IV in the original data is a factor
      # find the name of the IV in the orig data, which is assumed to be 2nd
      IVorig <- rownames(attributes(terms(model))$factors)[2] 
      
      
      if (is.factor(rawdata[,IVorig])) { # is the rawdata IV a factor?
        
        if (length(unique((rawdata[,IVorig]))) == 2) {
          #			if (length(unique((rawdata[,IV]))) == 2) {
          message('\nThe IV is a factor (not a continuous variable) with 2 levels.')
          message('\nThe values of 0 and 1 will be used for the 2 IV levels.')
          IV_min <- 0
          IV_max <- 1
        } else {message('\nThe IV is a factor (not a continuous variable) with more than 2 levels. 
				             Expect errors. The results are not meaningful.\n\n')
        }
      } else {
        IV_min <- min(rawdata[,IV])
        IV_max <- max(rawdata[,IV])
      }
    } else {
      IV_min <- IV_range[1]
      IV_max <- IV_range[2]
    }			
    
    
    # MOD values 
    if (is.null(MOD_range)) {
      if (is.factor(rawdata[,MOD])) { # is the rawdata MOD a factor?
        message('\nMOD is a factor, whereas it should be numeric. Expect errors.\n\n')
      }
      MOD_min <- min(rawdata[,MOD])
      MOD_max <- max(rawdata[,MOD])
      byint <- abs(MOD_min - MOD_max) / 100
      MODvalues <- seq(MOD_min, MOD_max, by = byint)
    } else {
      MOD_min <- MOD_range[1]
      MOD_max <- MOD_range[2]
      byint <- abs(MOD_min - MOD_max) / 100
      MODvalues <- seq(MOD_min, MOD_max, by = byint)
    }
    
    
    # get B & S
    if (is.null(names_IV_MOD)) {
      B <- fixef(model)              # B[1] = int,  B[2] = IV,  B[3] = MOD,  B[4] = product
      S <- stats::vcov(model)
    } else {
      B <- fixef(model)[c('(Intercept)',names_IV_MOD)]
      S <- stats::vcov(model)[c('(Intercept)',names_IV_MOD), c('(Intercept)',names_IV_MOD)]
    }
    
    # critical t with df of fixed effects 
    tcrit <- qt(0.975, model$fixDF$X[1])
    
    # construct the quadratic equation

    a <- tcrit^2 * S[PROD,PROD] - B[PROD]^2
    b <- 2 * (tcrit^2 * S[IV,PROD] - B[IV]*B[PROD])
    c <- tcrit^2 * S[IV,IV] - B[IV]^2
    
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
    message('\n   The ',IV,' values range from ',round(IV_min,2),' to ',round(IV_max,2), '\n\n')
  }
  
  

  
  if (!is.null(JN.data)) {
    
    if (is.null(Xaxis_label))  Xaxis_label <- MOD
    
    if (is.null(Yaxis_label))  Yaxis_label <- c(paste("Simple Slopes of",IV,'on',DV))
    
    if (is.null(plot_title))   plot_title <- c(paste("Simple Slopes of",IV,"on",DV,"by",MOD))
    
    if (is.null(legend_label))  legend_label <- 'Simple Slopes'
    
    ylimMIN <- min(JN.data[,1:3])
    ylimMAX <- max(JN.data[,1:3])
    
    MOD_min <- min(JN.data$MODvalues)
    MOD_max <- max(JN.data$MODvalues)
    
    plot(JN.data$MODvalues, JN.data$SimpleSlope, type='l', lty=1, col = "red", lwd = 2, cex.lab=1.3, 
         xlab=Xaxis_label, 
         ylab=Yaxis_label, 
         main=plot_title, 
         ylim=c(ylimMIN,ylimMAX) )
    
    # placing grey rectangle for region of NON sig	   
    if (ros[1] < MOD_min) Xboxmin <- MOD_min - abs(MOD_min) # to make sure the box starts at the y axis, if so
    if (ros[1] > MOD_min) Xboxmin <- ros[1]
    if (ros[2] > MOD_max) Xboxmax <- MOD_max + abs(MOD_max) # to make sure the box ends at the left of plot, if so
    if (ros[2] < MOD_max) Xboxmax <- ros[2]
    
    rect(Xboxmin, ylimMIN*3, Xboxmax, ylimMAX*3, col = 'gray93')
    
    lines(JN.data$MODvalues,  JN.data$SimpleSlope, type='l', lty=1, col = "red",  lwd = 2)
    lines(JN.data$MODvalues,  JN.data$CI.L,        type='l', lty=3, col = "red",  lwd = 1)
    lines(JN.data$MODvalues,  JN.data$CI.U,        type='l', lty=3, col = "red",  lwd = 1)
    
    # placing blue line for the lower region of sig	   
    if (ros[1] > MOD_min) abline(v=ros[1], col='blue', lty="solid",  lwd = 2)
    # placing blue line for the upper region of sig	   
    if (ros[2] < MOD_max) abline(v=ros[2], col='blue', lty="solid",  lwd = 2)
    
    abline(h=0, untf=FALSE,lty=3,lwd=1) 
    
    if (ros[1] < MOD_min & ros[2] < MOD_max) {
      legJNmax <- c(paste('J-N high >=',round(ros[2],2)))
      legend("top", legend = c(legend_label,'95% CI',legJNmax), 
             bty="n", lty=c(1,3,1), lwd=2, col=c('red','red','blue') )
    }	
    if (ros[1] > MOD_min & ros[2] > MOD_max) {
      legJNmin <- c(paste('J-N low <=',round(ros[1],2)))
      legend("top", legend = c(legend_label,'95% CI',legJNmin), 
             bty="n", lty=c(1,3,1), lwd=2, col=c('red','red','blue') )
    }
    if (ros[1] > MOD_min & ros[2] < MOD_max) {
      legJNmin <- c(paste('J-N low <=',round(ros[1],2)))
      legJNmax <- c(paste('J-N high >=',round(ros[2],2)))	
      legend("top", legend = c(legend_label,'95% CI',legJNmin,legJNmax), 
             bty="n", lty=c(1,3,1,1), lwd=2, col=c('red','red','blue','blue') )
    }
    
    ROSoutput <- list(JN.data = JN.data, ros=ros)
    class(ROSoutput) <- "MODERATED_REGRESSION"
    
    return(invisible(ROSoutput))
  }	
}
