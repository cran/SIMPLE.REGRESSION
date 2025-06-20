

SET_CORRELATION <- function(data, IVs, DVs, IV_covars=NULL, DV_covars=NULL,
                            Ncases=NULL, verbose=TRUE, display_cormats=FALSE) {
  
  if (anyDuplicated(DVs))        stop("There are duplicated elements in DVs")
  if (anyDuplicated(IVs))        stop("There are duplicated elements in IVs")
  if (anyDuplicated(DV_covars))  stop("There are duplicated elements in DV_covars")
  if (anyDuplicated(IV_covars))  stop("There are duplicated elements in IV_covars")
  if (anyDuplicated(c(DVs,IVs))) stop("There are duplicated elements in DVs & IVs")
  
  if (any(DVs %in% IVs))        stop("Some elements of DVs also appear in IVs")
  if (any(DVs %in% DV_covars))  stop("Some elements of DVs also appear in DV_covars")
  if (any(IVs %in% IV_covars))  stop("Some elements of IVs also appear in IV_covars")
  if (any(DVs %in% IV_covars))  stop("Some elements of DVs also appear in IV_covars")
  if (any(IVs %in% DV_covars))  stop("Some elements of IVs also appear in DV_covars")
  
  # is data a correlation matrix 
  if (nrow(data) == ncol(data)) {
    if (all(diag(data==1))) {datakind = 'cormat'}} else{datakind = 'raw'
    }
  
  if (datakind == 'cormat' & is.null(Ncases))  
    message('\nNcases must be provided when data is a correlation matrix.\n')
  
  if (datakind == 'cormat')  bigR <- data[x = c(IVs, IV_covars, DVs, DV_covars),
                                          y = c(IVs, IV_covars, DVs, DV_covars)]
  
  if (datakind == 'raw') {
    
    donnes <- data
    
    if (anyNA(donnes)) {
      
      donnes <- na.omit(donnes)
      
      cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
    }	
    
    Ncases <- nrow(donnes)
    
    bigR <- cor(donnes[,c(IVs, IV_covars, DVs, DV_covars)])
  }
  
  
  # use same names as in Cohen 1982 Table 1 p 315
  setD <- DVs
  setC <- DV_covars
  setB <- IVs
  setA <- IV_covars
  
  
  # covariance matrices from Cohen 1982 Table 1 p 315
  
  # Whole
  if (is.null(DV_covars) & is.null(IV_covars)) {
    
    type <- 'Whole'
    
    Ryx <- bigR[x = c(IVs, DVs), y = c(IVs, DVs)]
  }
  
  # Partial
  if (!is.null(DV_covars) & !is.null(IV_covars)) {
    if (all(DV_covars == IV_covars)) {
      
      type <- 'Partial'
      
      Rdb <- bigR[x = setD, y = setB]
      
      Rd  <- bigR[x = setD, y = setD]
      Ra  <- bigR[x = setA, y = setA]
      
      Rda <- bigR[x = setD, y = setA]
      
      Cd.a <- Rd - Rda %*% solve(Ra) %*% t(Rda) 
      
      Rb  <- bigR[x = setB, y = setB]
      Rba <- bigR[x = setB, y = setA]
      
      Cb.a <- Rb - Rba %*% solve(Ra) %*% t(Rba) 
      
      Rab <-  bigR[x = setA, y = setB]
      
      Cd.a_b.a <- Rdb - Rda %*% solve(Ra) %*% Rab 
      
      # cormat for analysis (IVs & DVs with covars partialled out)
      Ryx <- cov2cor(rbind( cbind(Cb.a, t(Cd.a_b.a)), cbind(Cd.a_b.a, Cd.a) ))
    }
  }
  
  # X Semipartial
  if (is.null(DV_covars) & !is.null(IV_covars)) {
    
    # note: the betas from this function, SYSTAT, and lmCor (psych) for X semipartial set correl
    # are all the same. However, the F tests from SYSTAT (for the 
    # individual DV regressions) are apparently wrong =
    # they are the same as the values for a full partial set correl from SYSTAT &  
    # the others. The F tests from this function & from lmCor are the 
    # same and are apparently the correct ones
    
    type <- 'X Semipartial'
    
    Rd <- bigR[x = setD, y = setD]
    
    Rb  <- bigR[x = setB, y = setB]
    Ra  <- bigR[x = setA, y = setA]
    Rba <- bigR[x = setB, y = setA]
    Rab <- bigR[x = setA, y = setB]
    Rdb <- bigR[x = setD, y = setB]
    
    Rda <- bigR[x = setD, y = setA]
    
    Cb.a <- Rb - Rba %*% solve(Ra) %*% t(Rba) 
    
    Cd_b.a <- Rdb - Rda %*% solve(Ra) %*% Rab 
    
    Ryx <- cov2cor(rbind( cbind(Cb.a, t(Cd_b.a)), cbind(Cd_b.a, Rd) ))
    
    H <- Cd_b.a %*% solve(Cb.a) %*% t(Cd_b.a)
    
    Rdab <- bigR[x = setD, y = c(setA, setB)]
    
    RAB <- bigR[x = c(setA,setB), y = c(setA,setB)]
    
    E <- Rd - Rdab %*% solve((RAB)) %*% t(Rdab)
    
    L <- det(E) / (det(E + H))
  }
  
  # Y Semipartial
  if (!is.null(DV_covars) & is.null(IV_covars)) {
    
    type <- 'Y Semipartial'
    
    Rd  <- bigR[x = setD, y = setD]
    Rc  <- bigR[x = setC, y = setC]
    Rdc <- bigR[x = setD, y = setC]
    
    Cd.c <- Rd - Rdc %*% solve(Rc) %*% t(Rdc)
    
    Rb <- bigR[x = setB, y = setB]
    
    Rdb <- bigR[x = setD, y = setB]
    
    Rcb <- bigR[x = setC, y = setB]
    
    Cd.c_b <- Rdb - Rdc %*% solve(Rc) %*% Rcb
    
    Ryx <- cov2cor(rbind( cbind(Rb, t(Cd.c_b)), cbind(Cd.c_b, Cd.c) ))
    
    Rbc <- bigR[x = setB, y = setC]
    
    Cb.c <- Rb - Rbc %*% solve(Rc) %*% t(Rbc) 
    
    # 1993 Cohen p 175 footnote: Says there was an error in the formula for H
    # in 1982 Cohen Table 2 for Y semipartial. But the Cd.c_b.c element in the 
    # formula he provided in 1993 still does not seem right. 
    # Using Cd.c_b instead works i.e., it reproduces the Rao F values from SYSTAT
    H <- Cd.c_b %*% solve(Cb.c) %*% t(Cd.c_b)
    
    E <- Cd.c - H
  }
  
  # Bipartial
  if (!is.null(DV_covars) & !is.null(IV_covars)) {
    if (all(DV_covars != IV_covars)) {
      
      type <- 'Bipartial'
      
      Rd  <- bigR[x = setD, y = setD]
      Rc  <- bigR[x = setC, y = setC]
      Rdc <- bigR[x = setD, y = setC]
      
      Cd.c <- Rd - Rdc %*% solve(Rc) %*% t(Rdc)
      
      Rb  <- bigR[x = setB, y = setB]
      Ra  <- bigR[x = setA, y = setA]
      Rba <- bigR[x = setB, y = setA]
      
      Cb.a <- Rb - Rba %*% solve(Ra) %*% t(Rba)
      
      Rdb <- bigR[x = setD, y = setB]
      Rcb <- bigR[x = setC, y = setB]
      Rcb <- bigR[x = setB, y = setC] # this version of Rcb works
      
      Rda <- bigR[x = setD, y = setA]
      Rab <- bigR[x = setA, y = setB]
      Rca <- bigR[x = setC, y = setA]
      Rca <- bigR[x = setA, y = setC] # this version of Rca works
      
      Cd.c_b.a <- Rdb - Rdc %*% solve(Rc) %*% t(Rcb) - Rda %*% solve(Ra) %*% Rab +
        Rdc %*% solve(Rc) %*% t(Rca) %*% solve(Ra) %*% Rab
      
      Ryx <- cov2cor(rbind( cbind(Cb.a, t(Cd.c_b.a)), cbind((Cd.c_b.a), Cd.c) ))
      
      H <- Cd.c_b.a %*% solve(Cb.a) %*% t(Cd.c_b.a)
      
      Rdab <- bigR[x = setD, y = c(setA, setB)]
      Rcab <- bigR[x = setC, y = c(setA, setB)]
      
      Cdc.ab <- Rdab - Rdc %*% solve(Rc) %*% Rcab
      
      RAB <- bigR[x = c(setA, setB), y = c(setA, setB)]
      
      E <- Cd.c - (Cdc.ab) %*% solve((RAB)) %*% t(Cdc.ab)
      
      L <- det(E) / (det(E + H))
    }
  }
  
  
  kx <- length(setB)
  ky <- length(setD)
  ka <- length(setA)
  kc <- length(setC)
  
  Rxx <- Ryx[1:kx,1:kx] 
  
  Ryy <- Ryx[(kx+1):ncol(Ryx),(kx+1):ncol(Ryx)] 
  
  Rx_y <- Ryx[1:kx,(kx+1):ncol(Ryx)]
  
  
  R_square <- 1 - ( det(Ryx) / (det(Ryy) * det(Rxx)) )
  
  
  # Rao's F
  
  kg <- 0  
  
  m <- Ncases - max(kc, ka+kg) - (ky + kx + 3) / 2
  
  s <- sqrt( (ky**2 * kx**2 - 4) / (ky**2 + kx**2 - 5))
  
  u <- ky * kx
  
  v <- m * s + 1 -u/2
  
  # Cohen 1988 p 470: "When Model I error is used, for all but the bipartial 
  # type of association, it can be shown that  L <- 1 - R_square
  # for Bipartial, use the formulas for E & H in Cohen 1982 p 320 Table 2, and
  # formula 10.1.3 from Cohen 1988 p 470
  if (type == 'Whole' | type == 'Partial') { 
    L <- 1 - R_square } else { L <- det(E) / (det(E + H)) }
  
  RaoF <- (L**(-1/s) - 1) * (v/u)
  
  
  # shrunken R_square -- Cohen 2003 p 616
  R_square_pop <- 1 - (1 - R_square) * ((v + u) / v)**s
  
  # 1984 Cohen, Nee  p 908
  Myx  <- solve(Ryy) %*% t(Rx_y) %*% solve(Rxx) %*% Rx_y
  q <- min(ky, kx)   # 1988 p 468
  
  # T-square    2003 p 611    1988 p 473 & 478 & 480
  T_square <- sum(diag(Myx)) / q
  # T_square <- L**(-1/s) - 1
  
  # shrunken T_square -- Cohen 2003 p 616   1984 Cohen, Nee  p 912
  w <- q * (Ncases - ky - kx - max(ka, kc) - 1)
  # issue with this part: ((w + u) / w)
  # Cohen 2003 p 616 & Cohen 1993 p  176  have it as ((w + u) / u), but this does not work right
  # 1984 Cohen, Nee  p 912  have it as ((w + u) / w) which works correctly
  T_square_pop <- 1 - (1 - T_square) * ((w + u) / w)
  
  # P-square    1993 p 170   2003 p 612
  P_square <- sum(diag(Myx)) / ky
  
  # shrunken P_square    1993 p 176    Cohen 2003 p 616   
  P_square_pop <- T_square_pop * (kx / ky)
  
  # # 2003 CCAW p 611   
  # # In the previous example
  # q = 3; Ncases = 97; ky = 5; kx = 3; ka = 4; kc = 4; u = 15; T_square = .110
  # w <- q * (Ncases - ky - kx - max(ka, kc) - 1)   # 252
  # w <- 3 * (97 - 5 - 3 - max(4,4) - 1)   # 252
  # T_square_pop <- 1 - (1 - T_square) * ((w + u) / w)   # 0.05702381
  # T_square_pop <- 1 - (1 - .110) * ((252 + 15) / 252); T_square_pop   # 0.05702381
  # P_square_pop <- T_square_pop * (kx / ky); P_square_pop   # 0.03421429
  # 
  # # Had there been six diagnostic groups instead of four
  # u = 25; kx = 6 - 1
  # q <- min(ky, kx)
  # w <- q * (Ncases - ky - kx - max(ka, kc) - 1)   # 410
  # T_square_pop <- 1 - (1 - T_square) * ((w + u) / w); T_square_pop   # 0.05573171
  # P_square_pop <- T_square_pop * (kx / ky); P_square_pop   
  # I get this 0.03343902 only for kx = 3, which is wrong
  
  
  # individual regressions (for each DV)
  
  # the commands below always replicate the R2 & betas from SYSTAT, but not always the F values
  # i.e., not for bipartial, Y semipartial, or X semipartial
  # (see the above note under Partial)
  # I tested the accuracy of the present F values in 2 ways:
  # 1 = I created residualized DVs & ran the lm function, which reproduced my R2 & F values
  # 2 = I used the residualized DVs in psych::lmCor, which also reproduced my R2 & F values
  
  # betas
  betas <- solve(Rxx) %*% Rx_y
  
  # Rsquared for each regression
  # R2 <- diag(t(Rx_y)  %*% solve(Rxx) %*% Rx_y)     same as:
  R2 <- colSums(betas * Rx_y)   # 2003 CCAW p 83
  
  # F for each regression   1984 Cohen, Nee  p 909 (?)
  df <- Ncases - kx - ka - 1

  F <- (R2 * df) / ((1 - R2) * kx)

  pF <- 1 - -expm1(pf(F, kx, df, lower.tail=FALSE, log.p=TRUE))
  # Fp <- pf(F, df1=m, df2=df, lower.tail=FALSE)
  
  R2Fmat <- data.frame(R2=R2, F=F, p=pF);   round(R2Fmat,6)
  
  # partial & semi-partial correlations -- need the partials for the t tests
  R_partials <- R_semipartials <- matrix(NA, nrow(Rx_y), ncol(Rx_y))
  N_preds <- length(IVs)
  for (lupeDVs in 1:length(DVs)) {
    partials <- matrix(NA,N_preds,1)
    # the DV must be in column 1 of cormat
    cormat <- Ryx[ c(DVs[lupeDVs], IVs), c(DVs[lupeDVs], IVs)]
    Rinv <- solve(cormat)
    for (luper in 2:(N_preds+1)) {
      partials[(luper-1),1] <- Rinv[luper,1] / ( -1 * sqrt(Rinv[1,1] * Rinv[luper,luper]) )
    }
    R_partials[,lupeDVs] <- partials
    
    semipartials <- sqrt( ( partials^2 / (1 - partials^2) ) * (1 - R2[lupeDVs]) )
    semipartials <- semipartials * sign(partials)
    
    R_semipartials[,lupeDVs] <- semipartials
    
    rownames(R_partials) <- rownames(R_semipartials) <- rownames(betas)
    colnames(R_partials) <- colnames(R_semipartials) <- colnames(betas)
  }  

  # https://online.stat.psu.edu/stat505/book/export/html/668
  t <- R_partials * sqrt( df / (1 - R_partials^2) )
  
  # this one comes close, but does not replicate SPSS t values:
  # t <- betas * sqrt( df / (1 - R2) )  # 2003 CCAW p 89
  
  se_betas <- betas / t   # boc:  just used  t = betas / se
  
  pt <- 2*pt(-abs(t), df=df)
  
  
  if (verbose) {
    
    cat('\n\n\nSet Correlation Analysis:\n')
    
    cat('\nThe DVs are:\n')
    # print(matrix(DVs,length(DVs),1), row.names = FALSE)
    cat('    ', matrix(DVs,length(DVs),1), sep = "\n    ")
    
    if (is.null(DV_covars)) { cat('\nNo variables were partialled out of the DVs.\n')
    } else {
      cat('\n\nThe variables partialled out of the DVs are:\n')
      # print(DV_covars)
      cat('    ', matrix(DV_covars,length(DV_covars),1), sep = "\n    ")
    }
    
    cat('\n\nThe IVs are:\n')
    # print(IVs)
    cat('    ', matrix(IVs,length(IVs),1), sep = "\n    ")
    
    if (is.null(IV_covars)) { cat('\nNo variables were partialled out of the IVs.\n')
    } else {
      cat('\n\nThe variables partialled out of the IVs are:\n')
      # print(IV_covars)
      cat('    ', matrix(IV_covars,length(IV_covars),1), sep = "\n    ")
    }
    
    cat('\n\nThe type of association is: ', type, '\n\n')
    
    if (display_cormats) {
      
      cat('\nUnpartialled correlations between all variables:\n\n')
      print(round(bigR,2))
      
      cat('\nCorrelations between the DVs after partialling out DV_covars (if any):\n\n')
      print(round(Ryy,3), print.gap=2)
      
      cat('\nCorrelations between the IVs after partialling out IV_covars (if any):\n\n')
      print(round(Rxx,3), print.gap=2)
    }
    
    cat('\n\nMeasures of multivariate association between the two sets:\n')
    cat('\n   R-square = ', round(R_square,3), '    Population R-square = ', round(R_square_pop,3))
    cat('\n   T-square = ', round(T_square,3))   #, '    Population T-square = ', round(T_square_pop,3))
    cat('\n   P-square = ', round(P_square,3))   #, '    Population P-square = ', round(P_square_pop,3))
    
    cat('\n\n\nRao multivariate F test of the null hypothesis:\n')
    cat('\n   Rao F = ', round(RaoF,3), '    df:  u = ', u, '   v = ', round(v,3))
    
    cat('\n\n\nRsquare and F test for the prediction of each basic y variable (df = ', kx, ', ', df,')\n\n', sep='')
    print(round(R2Fmat,3), print.gap=4)
    
    cat('\n\nCorrelations between the IVs & DVs after partialling out IV_covars (if any) & DV_covars (if any):\n\n')
    print(round(Rx_y,3), print.gap=2)
    
    cat('\n\nBetas predicting the DVs (columns) from the IVs (rows):\n\n')
    print(round(betas,3))
    
    cat('\n\nStandard errors of betas:\n')
    print(round(se_betas,3))
    
    cat('\n\nt-Statistic for betas:\n')
    print(round(t,3))
    
    cat('\n\np-Value for betas:\n')
    print(round(pt,5))
    
  }
  
  output <- list(bigR=bigR, Ryy=Ryy, Rxx=Rxx, Rx_y=Rx_y,
                 betas=betas, se_betas=se_betas, t=t, pt=pt)
  
  return(invisible(output))
  
  
  # measures of multivariate association between two sets of variables
  # R_square is interpretable as the proportion of generalized variance
  # ("multivariance") in set Y accounted for by set X
  
  # T-square  the proportion of the total variance of the smaller set accounted for by the larger.
  # the proportionof the total canonical variance that the sets account for in each other.
  
}
