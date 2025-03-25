

PLOT_MODEL <- function(model, 
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
  
  
  # type of model 
  if (inherits(model,"OLS_REGRESSION"))            model_type = 'OLS'
  if (inherits(model,"MODERATED_REGRESSION"))      model_type = 'MODERATED'
  if (inherits(model,"LOGISTIC_REGRESSION"))       model_type = 'LOGISTIC'
  if (inherits(model,"COUNT_REGRESSION"))          model_type = 'COUNT'
  
  if (any(grepl( 'pscl::zeroinfl',  model$modelMAIN$call, fixed = TRUE)))  model_type = 'ZINFL'
  
  # if (any(class(model$modelMAIN) == 'lm'))                          model_type = 'OLS'
  # if (any(grepl( 'binomial', model$modelMAIN$call, fixed = TRUE)))  model_type = 'logistic'
  # if (any(grepl( 'poisson',  model$modelMAIN$call, fixed = TRUE)))  model_type = 'poisson'
  
  if (inherits(model,"MODERATED_REGRESSION"))   names(model)[5] <- 'modelMAIN'
  
  
  # # yhatR$se.fit is not available for 'zinfl_poisson' & 'zinfl_negbin' models 
  # # so can't do regular CIs, must bootstrap
  # bootstrap_flag <- FALSE
  # if (!bootstrap & model_type == 'ZINFL')  bootstrap <- bootstrap_flag <- TRUE
  
  
  # DV name
  if (model_type == 'ZINFL') {
    DV <- names(attr(model$modelMAIN$terms$full, "dataClasses"))[1]
  } else { DV <- colnames(model$modeldata)[1] }
  
  # IV names
  if (model_type == 'ZINFL') {
    IVnames <- attr(model$modelMAIN$terms$full, "term.labels")
  } else { IVnames <- attr(model$modelMAIN$terms, "term.labels") }
  
  # remove interaction terms i.e., that contain :
  IVnames <- IVnames[!grepl(':', IVnames)]
  
  # are any of the predictors factors?
  if (model_type == 'ZINFL') { 
    list_xlevels <- model$modelMAIN$levels
  } else { list_xlevels <- model$modelMAIN$xlevels }
  
  
  
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
  
  # if IV_focal_1 is numeric & IV_focal_1_values were not provided
  if (is.null(IV_focal_1_values)) {
    
    IV_focal_1_range <- range(model$modeldata[IV_focal_1])
    
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
    
    # if IV_focal_2 is numeric & IV_focal_2_values were not provided -- have to use just a few values
    if (is.null(IV_focal_2_values)) {
      
      IV_focal_2_mn <- mean(unlist(model$modeldata[IV_focal_2]))
      
      IV_focal_2_sd <- sd(unlist(model$modeldata[IV_focal_2]))
      
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
        
      } else { IVs_nonfocal_values_list[[IVs_nonfocal[lupe]]] <- mean(model$modeldata[,IVs_nonfocal[lupe]]) }
      
    }
    
    # if IVs_nonfocal_values were provided by user, substitute them into IVs_nonfocal_values_list
    # IVs_nonfocal_values = list( EDUC = 5)
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
  
  # making sure that the order of the variables for testdata is the same as for model$modelMAIN
  testdata <- testdata[IVnames]
  
  head(testdata)
  
  # getting the predicted values & CIs
  testdata <- cbind(testdata, predict_boc(modelMAIN=model$modelMAIN, modeldata=model$modeldata,
                                          newdata=testdata, CI_level = 95, bootstrap=bootstrap, 
                                          model_type=model_type, family = model$family))
  
  # are CIs available?
  if (anyNA(testdata)) { CIs <- FALSE } else {CIs <- TRUE}
  
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (plot_save == TRUE) {
    
    plot_title = title
    
    if (!model_type == 'ZINFL') {height=7; width=9}
    if (model_type  == 'ZINFL') {height=9; width=7}
    
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
    
    # if (!model_type == 'ZINFL')
    #   par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)
  }
  
  
  
  if  (!model_type == 'ZINFL') {
    
    par(mfrow=c(1,1), pty="m", mar=c(3,2,3,2) + 2.6)
    
    # renaming yhat
    found <- match(colnames(testdata), 'yhat', nomatch = 0)
    DV_predicted <- paste(DV, '_predicted', sep='')
    colnames(testdata)[colnames(testdata) %in% 'yhat'] <- DV_predicted
    
    head(testdata)
    
    # if (is.null(xlim))  xlim <- c(min(testdata[,IV_focal_1]), max(testdata[,IV_focal_1]))
    
    if (is.null(xlab))  xlab <- IV_focal_1
    
    if (is.null(ylim)) {
      
      if (model_type == 'OLS' | model_type == 'MODERATED')  
        ylim <- c(min(testdata$ci_lb), max(testdata$ci_ub))
      
      if (model_type == 'LOGISTIC')  ylim <- c(0, 1)
      
      if (model_type == 'COUNT')  {
        
        if (CIs) {ylim <- c(0, max(testdata$ci_ub))
        } else {ylim <- c(0, max(model$modeldata[DV]))} # max(testdata[,DV_predicted])
      }
    }
    
    if (is.null(ylab)) {
      
      if (model_type == 'OLS' | model_type == 'MODERATED')  ylab <- DV
      
      if (model_type == 'LOGISTIC')  ylab <- paste("Probability of ", DV)
      
      if (model_type == 'COUNT')     ylab <- DV
    }
    
    if (is.null(title)) {
      
      if (model_type == 'OLS' | model_type == 'MODERATED') title <- paste('OLS regression prediction of', DV)
      if (model_type == 'LOGISTIC')  title <- paste('Logistic regression prediction of', DV)
      if (model_type == 'COUNT')     title <- paste('Count regression prediction of', DV)
    }  
    
    
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
        
        polygon(c(rev(testdata[,IV_focal_1]), testdata[,IV_focal_1]), 
                c(rev(testdata$ci_ub), testdata$ci_lb), col = 'grey95', border = NA)
        if (CIs) {
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
        
        legend("topleft", legend=legvalues, title=IV_focal_2, col=cols, # 1:length(IV_focal_2_values), 
               bty="n", lty=1, lwd=2)
        
        grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
      }
    } 
    
    
    # if IV_focal_1 is a factor
    if (IV_focal_1 %in% names(list_xlevels)) {
      
      # colours for the IV_focal_1_values
      if (length(cols_user) >= length(IV_focal_1_values))  {
        cols <- cols_user[1:length(IV_focal_1_values)]
      } else {cols <- rainbow(length(IV_focal_1_values))}
      
      if (is.null(ylim)) {
        if (CIs) {ylim <- range( pretty(testdata$ci_lb), pretty(testdata$ci_ub))
        } else {ylim <- range( pretty(testdata[,DV_predicted]), pretty(testdata[,DV_predicted]))}
      }
      
      if (!model_type == 'LOGISTIC')  ylim[2] <- ylim[2] + (ylim[2] * .05)
      
      if (is.null(IV_focal_2)) {  
        
        plot_bar <- barplot(t(testdata[DV_predicted]), beside=TRUE, legend.text=FALSE, col=cols,
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
        
        graphics::legend("top", legend = IV_focal_1_values, bty="n", lwd=4, col=cols, cex = .80)
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
  
  
  
  if  (model_type == 'ZINFL') {
    
    par(mfrow = c(2, 1),  mar = c(5,6,4,12))#, xpd= NA )  # oma = c(0,4,0,4) )#, mar = c(2,2,1,1)).  , xpd=TRUE
    
    if (is.null(xlab))  xlab <- IV_focal_1
    
    testdata_ZINFL <- testdata
    
    for (lupe12 in 1:2) {
      
      # zero data
      if (lupe12 == 1) {
        # remove the count data
        testdata <- testdata_ZINFL[,!grepl("count", colnames(testdata_ZINFL))]
        # remove "_zero" from column names
        names(testdata) <- gsub(pattern = "_zero", replacement = "", x = names(testdata))
        model_type = 'LOGISTIC'
        ylim <- ylab <- title <- NULL
      }
      
      # count data
      if (lupe12 == 2) {
        # remove the zero data
        testdata <- testdata_ZINFL[,!grepl("zero", colnames(testdata_ZINFL))]
        # remove "_count" from column names
        names(testdata) = gsub(pattern = "_count", replacement = "", x = names(testdata))
        model_type = 'COUNT'
        ylim <- ylab <- title <- NULL
      }
      
      # renaming yhat
      found <- match(colnames(testdata), 'yhat', nomatch = 0)
      DV_predicted <- paste(DV, '_predicted', sep='')
      colnames(testdata)[colnames(testdata) %in% 'yhat'] <- DV_predicted
      
      head(testdata)
      
      
      if (is.null(ylim)) {
        
        if (model_type == 'LOGISTIC')  ylim <- c(0,1)
        
        if (model_type == 'COUNT' | model_type == 'ZINFL')  {
          
          if (CIs) {ylim <- c(0, max(testdata$ci_ub))} else {ylim <- c(0, max(model$modeldata[DV]))} # max(testdata[,DV_predicted])
        }
      }
      
      if (is.null(ylab)) {
        
        if (model_type == 'LOGISTIC')  ylab <- paste("Probability of ", DV)
        
        if (model_type == 'COUNT' | model_type == 'ZINFL')   ylab <- DV
      }
      
      if (is.null(title)) {
        
        if (model_type == 'LOGISTIC')  title <- paste('Logistic regression prediction of', DV)
        
        if (model_type == 'COUNT')     title <- paste('Count regression prediction of', DV)
      }  
      
      
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
          
          grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
          
          
          for (lupe in 2:length(IV_focal_2_values)) {
            
            dum <- testdata[testdata[IV_focal_2] == IV_focal_2_values[lupe],]
            
            lines(dum[,IV_focal_1], dum[,DV_predicted], col = cols[lupe], lwd = 1.3)
          }
          
          legvalues <- IV_focal_2_values
          if (is.numeric(legvalues))  legvalues <- round(legvalues,3)
          
          legend("topright", legend=legvalues, title=IV_focal_2, col=cols, 
                 bty="n", lty=1, lwd=2, inset=c(-.4, .2), xpd= NA) 
          
          grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
        }
        
      } 
      
      
      # if IV_focal_1 is a factor
      if (IV_focal_1 %in% names(list_xlevels)) {
        
        # colours for the IV_focal_1_values
        if (length(cols_user) >= length(IV_focal_1_values))  {
          cols <- cols_user[1:length(IV_focal_1_values)]
        } else {cols <- rainbow(length(IV_focal_1_values))}
        
        if (is.null(ylim)) {
          if (CIs) {ylim <- range( pretty(testdata$ci_lb), pretty(testdata$ci_ub))
          } else {ylim <- range( pretty(testdata$yhat), pretty(testdata$yhat))}
        }
        
        if (!model_type == 'LOGISTIC')  ylim[2] <- ylim[2] + (ylim[2] * .05)
        
        if (is.null(IV_focal_2)) {  
          
          plot_bar <- barplot(t(testdata[DV_predicted]), beside=TRUE, legend.text=FALSE, col=cols,
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
          
          graphics::legend("top", legend = IV_focal_1_values, bty="n", lwd=4, col=cols, cex = .80)
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
    
    testdata <- testdata_ZINFL
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

