

PLOT_MODEL <- function(modobject, 
                       IV_focal_1, IV_focal_1_values=NULL, 
                       IV_focal_2=NULL, IV_focal_2_values=NULL, 
                       IVs_nonfocal_values = NULL,
                       bootstrap=FALSE, N_sims=100, CI_level=95, 
                       xlim=NULL, xlab=NULL,
                       ylim=NULL, ylab=NULL,
                       title = NULL,
                       plot_save = FALSE, plot_save_type = 'png',
                       cols_user = NULL,
                       verbose=TRUE) {
  
  
  if (is.null(cols_user))
    cols_user <- c("blue", "red", 'black', 'cyan2', "blueviolet", 'limegreen', "yellow", 'mediumvioletred')
  # cols_user <- c("mediumvioletred", 'black', "blue", 'cyan2', "red", 'limegreen', "yellow", 'blueviolet')
  # cols_user <- c("ivory", 'black', "blue", 'cyan2', "red", 'limegreen', "yellow", 'blueviolet')
  
  
  # kind of modobject 
  if (inherits(modobject,"OLS_REGRESSION"))            kind = 'OLS'
  if (inherits(modobject,"MODERATED_REGRESSION"))      kind = 'MODERATED'
  if (inherits(modobject,"LOGISTIC_REGRESSION"))       kind = 'LOGISTIC'
  if (inherits(modobject,"COUNT_REGRESSION"))          kind = modobject$kind
  
  # model variable names
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    modvars <- names(attr(modobject$model$terms$full, "dataClasses"))
  } else { modvars <- names(attr(modobject$model$terms, "dataClasses")) }
  # modvars <- do.call("c", modobject$noms_list)
  
  DV <- modvars[1]
  
  IVnames <- modvars[-1]
  # # remove offset variable (count regression) if there was one
  # IVnames <- IVnames[!grepl('offset', IVnames)]
  IVnames <- gsub("offset\\(", "", IVnames)
  IVnames <- gsub("\\)", "", IVnames)
  
  modeldata <- modobject$modeldata[,c(DV, IVnames)]
  
  # remove interaction terms i.e., that contain :
  IVnames <- IVnames[!grepl(':', IVnames)]
  
  # are any of the predictors factors?
  if (kind == 'ZINFL' | kind == 'HURDLE') { 
    list_xlevels <- modobject$model$levels
  } else { list_xlevels <- modobject$model$xlevels }
  
  
  
  # IV_focal_1_values
  
  # test if IV_focal_1 is in the IVnames
  if (!is.element(IV_focal_1, IVnames)) 
    message('\nThe IV_focal_1 variable does not appear in the model predictors. Perhaps there is a spelling error.\n')
  
  # if IV_focal_1 is a factor
  if (IV_focal_1 %in% names(list_xlevels)) {
    
    # if IV_focal_1 is a factor & IV_focal_1_values were NOT provided, use all levels of it
    if (is.null(IV_focal_1_values))  IV_focal_1_values <- unlist(list_xlevels[IV_focal_1])
    
    # if IV_focal_1 is a factor & IV_focal_1_values were provided, test if valid levels of it were provided
    if (!is.null(IV_focal_1_values) & !all(is.element(IV_focal_1_values,  unlist(list_xlevels[IV_focal_1])))) {
      
      eles_false <- which(is.element(IV_focal_1_values, unlist(list_xlevels[IV_focal_1])) == FALSE)
      
      invalid_levels <- IV_focal_1_values[eles_false]
      
      message('\nThe following specified levels of "', IV_focal_1 , '" do not exist in the model: ', invalid_levels)
      
      IV_focal_1_values <- IV_focal_1_values [! IV_focal_1_values %in% invalid_levels]
    }
  }
  
  # if IV_focal_1 is not a factor (so must be numeric) & IV_focal_1_values were not provided
  if (!(IV_focal_1 %in% names(list_xlevels)) & is.null(IV_focal_1_values)) {
    
    IV_focal_1_range <- range(modeldata[IV_focal_1])
    
    IV_focal_1_values <- seq(IV_focal_1_range[1], IV_focal_1_range[2], length.out = 100)
  }
  
  
  
  # IV_focal_2_values
  
  # notice if IV_focal_2 is not provided but IV_focal_2_values are provided
  if (is.null(IV_focal_2) & !is.null(IV_focal_2_values)) 
    message('\nIV_focal_2_values were provided but without specifying the name of IV_focal_2.')
  
  if (!is.null(IV_focal_2)) {
    
    # test if IV_focal_2 is in the IVnames
    if (!is.element(IV_focal_2, IVnames)) 
      message('\nThe IV_focal_2 variable does not appear in the model predictors. Perhaps there is a spelling error.')
    
    # if IV_focal_2 is a factor
    if (IV_focal_2 %in% names(list_xlevels)) {
      
      # if IV_focal_2 is a factor & IV_focal_2_values were NOT provided, use all levels of it
      if (is.null(IV_focal_2_values))  IV_focal_2_values <- unlist(list_xlevels[IV_focal_2])
      
      # if IV_focal_2 is a factor & IV_focal_2_values were provided, test if valid levels of it were provided
      if (!is.null(IV_focal_2_values) & !all(is.element(IV_focal_2_values,  unlist(list_xlevels[IV_focal_2])))) {
        
        eles_false <- which(is.element(IV_focal_2_values, unlist(list_xlevels[IV_focal_2])) == FALSE)
        
        invalid_levels <- IV_focal_2_values[eles_false]
        
        message('\nThe following specified levels of "', IV_focal_2 , '" do not exist in the model: ', invalid_levels)
        
        IV_focal_2_values <- IV_focal_2_values [! IV_focal_2_values %in% invalid_levels]
      }
    }
    
    # if IV_focal_2 is not a factor (so must be numeric)  & IV_focal_2_values 
    # were not provided -- have to use just a few values
    if (!(IV_focal_2 %in% names(list_xlevels)) & is.null(IV_focal_2_values)) {
      
      IV_focal_2_mn <- mean(unlist(modeldata[IV_focal_2]))
      
      IV_focal_2_sd <- sd(unlist(modeldata[IV_focal_2]))
      
      IV_focal_2_values <- c((IV_focal_2_mn - IV_focal_2_sd), IV_focal_2_mn, (IV_focal_2_mn + IV_focal_2_sd))
    }
  }
  
  
  # non focal IV values
  
  IVs_nonfocal <- setdiff(IVnames, c(IV_focal_1, IV_focal_2))
  
  # remove interaction terms i.e., that contain :
  IVs_nonfocal <- IVs_nonfocal[!grepl(':', IVs_nonfocal)]
  
  IVs_nonfocal_values_list <- c()
  
  # when there are factors, cycle through the IVs_nonfocal, use mean if continuous, use baseline if categorical
  if (length(IVs_nonfocal) > 0) {
    for (lupe in 1:length(IVs_nonfocal)) {
      
      if (IVs_nonfocal[lupe] %in% names(list_xlevels)) {
        
        IVs_nonfocal_values_list[[IVs_nonfocal[lupe]]] <- sort(list_xlevels[[IVs_nonfocal[lupe]]])[1]
        
      } else { IVs_nonfocal_values_list[[IVs_nonfocal[lupe]]] <- mean(modeldata[,IVs_nonfocal[lupe]]) }
      
    }
    
    # if IVs_nonfocal_values were provided by user, substitute them into IVs_nonfocal_values_list
    if (!is.null(IVs_nonfocal_values)) {
      
      # first check to make sure the names of the IVs_nonfocal_values are in the data
      wrongnoms <- setdiff( names(IVs_nonfocal_values), names(IVs_nonfocal_values_list) )
      if (!identical(wrongnoms, character(0))) 
        message('\n\nOne or more names on IVs_nonfocal_values does not exist in the model data. Expect errors.\n')
      
      for (lupe in 1:length(IVs_nonfocal_values))
        IVs_nonfocal_values_list[names(IVs_nonfocal_values)[lupe]] <- IVs_nonfocal_values[lupe]
    }
  }
  
  IVs_nonfocal_values_df <- data.frame(cbind(lapply(IVs_nonfocal_values_list, 
                                                    function(y) if(is.numeric(y)) round(y, 2) else y)) )
  colnames(IVs_nonfocal_values_df) <- NULL
  
  
  
  # testdata for predict
  testdata <- IVs_nonfocal_values_list
  
  testdata[[IV_focal_1]] <- IV_focal_1_values
  
  if (!is.null(IV_focal_2)) testdata[[IV_focal_2]] <- IV_focal_2_values
  
  testdata <- do.call(expand.grid, testdata)
  
  # making sure that the order of the variables for testdata is the same as for modobject$model
  testdata <- testdata[IVnames]
  
  head(testdata)
  
  # getting the predicted values & CIs
  testdata <- cbind(testdata, predict_boc(model=modobject$model, modeldata=modeldata,
                                          newdata=testdata, CI_level=CI_level, bootstrap=bootstrap, 
                                          kind=kind, family=modobject$family))
  
  # are CIs available?
  # if (anyNA(testdata)) { CIs <- FALSE } else {CIs <- TRUE}
  CIs <- TRUE
  if (all(sapply(testdata, anyNA)) == FALSE)  CIs <- FALSE
  

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  # on.exit(par(...,, add = TRUE))
  
  
   
  
  if (plot_save == TRUE) {
    
    plot_title = title
    
    if (kind != 'ZINFL' & kind != 'HURDLE') {height=7; width=9}
    if (kind  == 'ZINFL' | kind == 'HURDLE')  {height=9; width=7}
    
    if (is.null(plot_save_type))  plot_save_type = 'png'
    
    if (plot_save_type == 'bitmap')
      bitmap(paste("Figure - ",plot_title,".bitmap",sep=""), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'tiff')
      tiff(paste("Figure - ",plot_title,".tiff",sep=""), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'png')
      png(paste("Figure - ",plot_title,".png",sep=""), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'jpeg')
      jpeg(paste("Figure - ",plot_title,".jpeg",sep=""), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'bmp')
      bmp(paste("Figure - ",plot_title,".bmp",sep=""), height=height, width=width, units='in', res=1200, pointsize=12)
    
    # if (kind != 'ZINFL' & kind != 'HURDLE')
    #   par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)
  }
  
  
  
  if (kind == 'OLS' | kind == 'MODERATED' | kind == 'LOGISTIC' | 
      kind == 'POISSON' | kind == 'NEGBIN') {
    
    par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)
    
    # renaming yhat
    found <- match(colnames(testdata), 'yhat', nomatch = 0)
    DV_predicted <- paste(DV, '_predicted', sep='')
    colnames(testdata)[colnames(testdata) %in% 'yhat'] <- DV_predicted
    
    head(testdata)
    
    if (is.null(xlim) & is.numeric((testdata[,IV_focal_1])))  
      xlim <- c(min(testdata[,IV_focal_1]), max(testdata[,IV_focal_1]))
    
    if (is.null(xlab))  xlab <- IV_focal_1
    
    if (is.null(ylim)) {
      
      if (kind == 'OLS' | kind == 'MODERATED')  
        # ylim <- c(min(testdata$ci_lb), max(testdata$ci_ub))
        ylim <- range( pretty(testdata$ci_lb), pretty(testdata$ci_ub))
      
      if (kind == 'LOGISTIC')  ylim <- c(0, 1)
      
      if (kind == 'POISSON' | kind == 'NEGBIN')  {
        
        # if (CIs) { ylim <- range(pretty(c(0, max(testdata$ci_ub))))   
        # } else {   ylim <- range(pretty(c(0, max(testdata[DV_predicted])))) }  
        
        if (CIs) { ylim <- c(0, max(testdata$ci_ub))
        } else {   ylim <- c(0, max(testdata[DV_predicted])) }
        
        # increase the upper ylim by 50%
        ylim <- better_ylim(ylim, buffer = 0.50)
        ylim[1] <- 0
      }
     }
    
    if (is.null(ylab)) {
      
      if (kind == 'OLS' | kind == 'MODERATED')  ylab <- DV
      
      if (kind == 'LOGISTIC')  ylab <- paste("Probability of ", DV)
      
      if (kind == 'POISSON' | kind == 'NEGBIN')     ylab <- DV
    }
    
    if (is.null(title)) {
      
      if (kind == 'OLS' | kind == 'MODERATED') title <- paste('OLS regression prediction of', DV)
      if (kind == 'LOGISTIC')  title <- paste('Logistic regression prediction of', DV)
      if (kind == 'POISSON' | kind == 'NEGBIN') title <- paste('Count regression prediction of', DV)
    }  
    
    plotfun(testdata=testdata, list_xlevels=list_xlevels, DV_predicted=DV_predicted, 
            CIs=CIs, kind=kind, cols_user=cols_user,
            xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
            IV_focal_1=IV_focal_1, IV_focal_1_values=IV_focal_1_values, 
            IV_focal_2=IV_focal_2, IV_focal_2_values=IV_focal_2_values)
  }
  
  
  if (kind == 'ZINFL' | kind == 'HURDLE') {
    
    # par(mfrow = c(2, 1),  mar = c(5,6,4,12))  #, xpd= NA )  # oma = c(0,4,0,4) )#, mar = c(2,2,1,1)).  , xpd=TRUE
    # par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)
    
    par(mfrow = c(2, 1),  oma = c(0,4,0,6) )  #, mar = c(2,2,1,1))
    
    testdata_ORIG <- testdata
    
    
    if (is.null(xlim) & !is.factor(testdata[,IV_focal_1]))  
      xlim <- c(min(testdata[,IV_focal_1]), max(testdata[,IV_focal_1]))
    
    if (is.factor(testdata[,IV_focal_1]))  xlim <- NULL
    
    if (is.null(xlab))  xlab <- IV_focal_1
    
    
    # zero data
    # remove the count data
    testdata <- testdata_ORIG[,!grepl("count", colnames(testdata_ORIG))]
    # remove "_zero" from column names
    names(testdata) <- gsub(pattern = "_zero", replacement = "", x = names(testdata))
    
    # reversing yhat for hurdle models, to make DV = the probability of crossing the hurdle
    if (kind == 'HURDLE') testdata$yhat <- 1 - testdata$yhat
    
    # renaming yhat
    found <- match(colnames(testdata), 'yhat', nomatch = 0)
    DV_predicted <- paste(DV, '_predicted', sep='')
    colnames(testdata)[colnames(testdata) %in% 'yhat'] <- DV_predicted
    
    ylim <- c(0, 1)
    
    ylab <- 'Probability'
    
    if (kind == 'ZINFL')  title <- paste(DV, ':\nProbability of excess zeroes', sep='')
    # if (kind == 'HURDLE') title <- paste(DV, ':\nProbability of zero', sep='')
    if (kind == 'HURDLE') title <- paste(DV, ':\nProb. of hurdle cross', sep='')
    
    
    plotfun(testdata=testdata, list_xlevels=list_xlevels, DV_predicted=DV_predicted, 
            CIs=CIs, kind=kind, cols_user=cols_user,
            xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
            IV_focal_1=IV_focal_1, IV_focal_1_values=IV_focal_1_values, 
            IV_focal_2=IV_focal_2, IV_focal_2_values=IV_focal_2_values)
    
    
    # count data
    # remove the zero data
    testdata <- testdata_ORIG[,!grepl("zero", colnames(testdata_ORIG))]
    # remove "_count" from column names
    names(testdata) = gsub(pattern = "_count", replacement = "", x = names(testdata))
    
    # renaming yhat
    found <- match(colnames(testdata), 'yhat', nomatch = 0)
    DV_predicted <- paste(DV, '_predicted', sep='')
    colnames(testdata)[colnames(testdata) %in% 'yhat'] <- DV_predicted
    
    if (CIs) { ylim <- c(0, max(testdata$ci_ub))
    } else {   ylim <- c(0, max(testdata[DV_predicted])) }
    
    # increase the upper ylim by 50%
    ylim <- better_ylim(ylim, buffer = 0.50)
    ylim[1] <- 0

    ylab <- DV
    
    title <- paste(DV, ':\nExpected counts', sep='')
    
    plotfun(testdata=testdata, list_xlevels=list_xlevels, DV_predicted=DV_predicted, 
            CIs=CIs, kind=kind, cols_user=cols_user,
            xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
            IV_focal_1=IV_focal_1, IV_focal_1_values=IV_focal_1_values, 
            IV_focal_2=IV_focal_2, IV_focal_2_values=IV_focal_2_values)
    
    testdata <- testdata_ORIG
  }
  
  
  if (plot_save == TRUE)  dev.off()
  
  
  if (verbose) {
    
    message('\nThe DV is "', DV, '"')
    
    # if IV_focal_1 is a factor vs numeric
    if (IV_focal_1 %in% names(list_xlevels)) { IV_focal_1_type <- 'a factor'
    } else { IV_focal_1_type <- 'numeric' }
    
    message('\nThe focal, varying IV is "', IV_focal_1, '", and it is ', IV_focal_1_type)
    
    if (!is.null(IV_focal_2)) {
      # if IV_focal_2 is a factor vs numeric
      if (IV_focal_2 %in% names(list_xlevels)) { IV_focal_2_type <- 'a factor'
      } else { IV_focal_2_type <- 'numeric' }
      
      message('\nThe second varying IV (moderator) is "', IV_focal_2, '," and it is ', IV_focal_2_type)
      message('\nThe levels of "', IV_focal_2, '" are:')
      if (is.numeric(IV_focal_2_values))  print(round(IV_focal_2_values,3))
      if (!is.numeric(IV_focal_2_values)) {
        rownames(IV_focal_2_values) <- NULL
        # print(unname(IV_focal_2_values), row.names = F)
        print(  unname(as.data.frame(IV_focal_2_values)), quote = FALSE, row.names = FALSE)
      }
    }
    
    # if (!is.null(IVs_nonfocal)) {
    if (length(IVs_nonfocal) > 0) {
      message('\nThe constant values for the other predictors are:')
      print(IVs_nonfocal_values_df, print.gap=4)
    }
    
    message('\nThe confidence interval percentage is: ', CI_level)
    
    message('\nBootstrapped confidence intervals: ', bootstrap)
    
    if (bootstrap) {
      # if (bootstrap_flag) {
      #   message('\nConventional CIs cannot be computed for zero-inflated models at this time.')
      #   message('Bootrapped CIs were computed instead.\n')
      # }
      message('\nThe number of bootstrap simulations: ', N_sims)
    }
    
    message('\nThe top rows of the output data matrix:\n')
    print(head(round_boc(testdata, 3), n=6), print.gap=4)
    
    message('\nThe bottom rows of the output data matrix:\n')
    print(tail(round_boc(testdata, 3), n=6), print.gap=4)
  }
  
  return(invisible(testdata)) 
}

