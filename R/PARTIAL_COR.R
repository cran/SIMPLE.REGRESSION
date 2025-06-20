

PARTIAL_COR <- function(data, Y, X=NULL, C=NULL, Ncases=NULL, verbose=TRUE) {
  
  if (as.character(match.call()[[1]]) == "PARTIAL_COEFS") {
    warning('\nA package function have been changed:
             \nplease use PARTIAL_COR instead of PARTIAL_COEFS for partial correlations.\n\n', call. = FALSE) 
  }
  
  if (anyDuplicated(Y))        stop("There are duplicated elements in Y")
  if (anyDuplicated(X))        stop("There are duplicated elements in X")
  if (anyDuplicated(C))        stop("There are duplicated elements in C")
  if (anyDuplicated(c(Y,X,C))) stop("There are duplicated elements across Y, X, & C")
  
  if (any(Y %in% X))  stop("Some elements of Y also appear in X")
  if (any(Y %in% C))  stop("Some elements of Y also appear in C")
  if (any(X %in% C))  stop("Some elements of X also appear in C")
  
  bigR <- Ry.c <- Rx_y <- betas <- R_partials <- R_semipartials <- se_betas <- t <- pt <- NA
  
  # is data a correlation matrix 
  if (nrow(data) == ncol(data)) {
    if (all(diag(data==1))) {datakind = 'cormat'}} else{datakind = 'raw'
    }
  
  if (datakind == 'cormat' & is.null(Ncases))  
    message('\nNcases must be provided when data is a correlation matrix.\n')         # FIX
  
  if (datakind == 'cormat')  bigR <- data[x = c(X, C, Y), y = c(X, C, Y)]
  
  if (datakind == 'raw') {
    
    donnes <- data[,c(X, C, Y)]
    
    if (anyNA(donnes)) {
      
      donnes <- na.omit(donnes)
      
      cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
    }	
    
    # check if all are numeric
    if (any(apply(donnes, 2, is.numeric) == FALSE)) 
      stop("Not all of the variabels are numeric.")
    
    Ncases <- nrow(donnes)
    
    bigR <- cor(donnes[,c(X, Y, C)])
  }
  
  if (!is.null(C) & length(Y) == 1) 
    cat('\n\nThere is only one Y variable. There are no Y correlations to be partialled by C.\n\n')
  
  if (!is.null(C) & length(Y) > 1) {
    
    Ryy <- bigR[x = Y, y = Y]
    
    if (length(C) > 1) { Rcc  <- bigR[x = C, y = C] } else { Rcc <- 1 }
    
    Ryc <- bigR[x = Y, y = C]
    
    Cy.c <- Ryy - Ryc %*% solve(Rcc) %*% t(Ryc) 
    
    Ry.c <- cov2cor(Cy.c)
    
    
    # # sample number
    # n <- dim(x)[1]
    # 
    # # given variables' number
    # gp <- dim(x)[2]-2
    # 
    # statistic <- pcor*sqrt((n-2-gp)/(1-pcor^2))
    # p.value <- 2*pt(-abs(statistic),(n-2-gp))
    # #p.value <- 2*pnorm(-abs(statistic))
    
    df <- Ncases - length(C) - 2
    
    t <- Ry.c * sqrt( df / (1 - Ry.c^2))
    
    pt <- 2*pt(-abs(t), df=df)
  }
  
  
  if (!is.null(X)) {
    
    # betas
    Rxx <- bigR[x = X, y = X] 
    
    Rx_y <- bigR[x = X, y = Y, drop=FALSE]
    
    betas <- solve(Rxx) %*% Rx_y
    
    # Rsquared for each regression
    # R2 <- diag(t(Rx_y)  %*% solve(Rxx) %*% Rx_y)     same as:
    R2 <- colSums(betas * Rx_y)   # 2003 CCAW p 83
    
    # partial & semi-partial correlations
    R_partials <- R_semipartials <- matrix(NA, nrow(Rx_y), ncol(Rx_y))
    N_preds <- length(X)
    for (lupeY in 1:length(Y)) {
      partials <- matrix(NA,N_preds,1)
      # the DV must be in column 1 of cormat
      cormat <- bigR[ c(Y[lupeY], X), c(Y[lupeY], X)]
      Rinv <- solve(cormat)
      for (luper in 2:(N_preds+1)) {
        partials[(luper-1),1] <- Rinv[luper,1] / ( -1 * sqrt(Rinv[1,1] * Rinv[luper,luper]) )
      }
      R_partials[,lupeY] <- partials
      
      semipartials <- sqrt( ( partials^2 / (1 - partials^2) ) * (1 - R2[lupeY]) )
      semipartials <- semipartials * sign(partials)
      
      R_semipartials[,lupeY] <- semipartials
      
      rownames(R_partials) <- rownames(R_semipartials) <- rownames(betas)
      colnames(R_partials) <- colnames(R_semipartials) <- colnames(betas)
      
      df <- Ncases - N_preds - 1
      
      # https://online.stat.psu.edu/stat505/book/export/html/668
      t <- R_partials * sqrt( df / (1 - R_partials^2))
      
      pt <- 2*pt(-abs(t), df=df)
      
      # t <- betas * sqrt( df / (1 - R2) )  # 2003 CCAW p 89
      # pt <- 2*pt(-abs(t), df=df)
      
      se_betas <- betas / t   # boc:  just used  t = betas / se
    }
  }
  
  if (verbose) {
    
    # cat('\nThe Y variables are:\n')
    # cat('    ', matrix(Y,length(Y),1), sep = "\n    ")
    
    if (!is.null(C) & length(Y) > 1) {
      # cat('\n\nThe C variables partialled out of the Y variables are:\n')
      # cat('    ', matrix(C,length(C),1), sep = "\n    ")
      
      cat('\n\nCorrelations between the Y & C variables:\n\n')
      print(round(Ryc,3), print.gap=2)
      
      cat('\n\nCorrelations between the Y variables before partialling:\n\n')
      print(round(Ryy,3), print.gap=2)
      
      cat('\n\nCorrelations between the Y variables after partialling out the C variables:\n\n')
      print(round(Ry.c,3), print.gap=2)
      
      cat('\n\nt-Statistic:\n')
      print(round(t,3))
      
      cat('\n\np-Values:\n')
      print(round(pt,5))
    }
    
    if (!is.null(X)) {
      
      # cat('\n\nThe X variables (the IVs) are:\n')
      # cat('    ', matrix(X,length(X),1), sep = "\n    ")
      
      if (length(X) == 1 & length(Y) == 1) {
        
        outmat <- rbind(Rx_y, betas, R_partials, R_semipartials, se_betas, t, pt)
        
        rownames(outmat) <- c('Rx_y','betas','R_partials','R_semipartials','se_betas','t','pt')
        
        print(round(outmat,3), print.gap=2)
      }
      
      if (length(X) > 1 & length(Y) == 1) {
        
        outmat <- cbind(Rx_y, betas, R_partials, R_semipartials, se_betas, t, pt)
        
        colnames(outmat) <- c('Rx_y','betas','R_partials','R_semipartials','se_betas','t','pt')
        
        print(round(outmat,3), print.gap=2)
      }
      
      if (length(X) > 1 & length(Y) > 1) {
        
        cat('\n\nBivariate correlations between the X & Y variables:\n\n')
        print(round(Rx_y,3), print.gap=2)
        
        cat('\n\nBetas for the X variables predicting the Y variables (columns):\n\n')
        print(round(betas,3), print.gap=2)
        
        cat('\n\nStandard errors of betas:\n')
        print(round(se_betas,3))
        
        cat('\n\nPartial correlations for the X variables predicting the Y variables:\n\n')
        print(round(R_partials,3), print.gap=2)
        
        cat('\n\nSemi-partial correlations for the X variables predicting the Y variables:\n\n')
        print(round(R_semipartials,3), print.gap=2)
        
        cat('\n\nt-Statistic for betas:\n')
        print(round(t,3))
        
        cat('\n\np-Value for betas:\n')
        print(round(pt,5))
      }
    }
  }
  
  output <- list(bigR=bigR, Ry.c=Ry.c, Rx_y=Rx_y, betas=betas, R_partials=R_partials,
                 R_semipartials=R_semipartials, se_betas=se_betas, t=t, pt=pt)
  
  return(invisible(output))
}

PARTIAL_COEFS <- PARTIAL_COR

