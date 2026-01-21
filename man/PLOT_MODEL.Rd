\name{PLOT_MODEL}
\alias{PLOT_MODEL}
\title{Plots predicted values for a regression model}
\description{Plots predicted values of the outcome variable for specified levels 
of predictor variables for OLS_REGRESSION, MODERATED_REGRESSION, 
LOGISTIC_REGRESSION, and COUNT_REGRESSION models from this package.}
\usage{PLOT_MODEL(modobject, 
           IV_focal_1, IV_focal_1_values=NULL, 
           IV_focal_2=NULL, IV_focal_2_values=NULL, 
           IVs_nonfocal_values = NULL,
           bootstrap=FALSE, N_sims=100, CI_level=95, 
           xlim=NULL, xlab=NULL,
           ylim=NULL, ylab=NULL,
           title = NULL,
           plot_save = FALSE, plot_save_type = 'png',
           cols_user = NULL,
           verbose=TRUE)
}
\arguments{
  \item{modobject}{
  \code{}The returned output from the OLS_REGRESSION, MODERATED_REGRESSION, 
  LOGISTIC_REGRESSION, or COUNT_REGRESSION functions in this package.}

  \item{IV_focal_1}{
  \code{}The name of the focal, varying predictor variable. 
  \cr Example: IV_focal_1 = 'age'}

  \item{IV_focal_1_values}{
  \code{}(optional) Values for IV_focal_1, for which predictions of the
  outcome will be produced and plotted.
  IV_focal_1_values will appear on the x-axis in the plot.
  If IV_focal_1 is numeric and IV_focal_1_values is not provided,
  then a sequence based on the range of the model data values for IV_focal_1 will be used.
  If IV_focal_1 is a factor & IV_focal_1_values is not provided, then the
  factor levels from the model data values for IV_focal_1 will be used.
  \cr Example: IV_focal_1_values = seq(20, 80, 1)
  \cr Example: IV_focal_1_values = c(20, 40, 60)}

  \item{IV_focal_2}{
  \code{}(optional) If desired, the name of a second focal predictor variable for the plot. 
  \cr Example: IV_focal_2 = 'height'}

  \item{IV_focal_2_values}{
  \code{}(optional) Values for IV_focal_2 for which predictions of the
  outcome will be produced and plotted.
  If IV_focal_2 is numeric and IV_focal_2_values is not provided, then
  the following three values for IV_focal_2_values, derived from the model data,
  will be used for plotting: the mean, one SD below the mean, and one SD above the mean.
  If IV_focal_2 is a factor & IV_focal_2_values is not provided, then the
  factor levels from the model data values for IV_focal_2 will be used.
  \cr Example: IV_focal_2_values = c(20, 40, 60)}

  \item{IVs_nonfocal_values}{
  \code{}(optional) A list with the desired constant values for the non focal predictors, 
  if any. If IVs_nonfocal_values is not provided, then the mean values of numeric non focal
  predictors and the baseline values of factors will be used as the defaults.
  It is also possible to specify values for only some of the IVs_nonfocal variables
  on this argument.
  \cr Example: IVs_nonfocal_values = list(AGE = 25, EDUC = 12)
 }

  \item{bootstrap}{
  \code{}(optional) Should bootstrapping be used for the confidence intervals?
  The options are TRUE or FALSE (the default).}

  \item{N_sims}{
  \code{}(optional) The number of bootstrap simulations. 
  \cr Example: N_sims = 1000}

  \item{CI_level}{
  \code{}(optional) The desired confidence interval, in whole numbers.
  \cr Example: CI_level = 95}

  \item{xlim}{
  \code{}(optional) The x-axis limits for the plot. 
  \cr Example: xlim = c(1, 9)}

  \item{xlab}{
  \code{}(optional) A x-axis label for the plot. 
  \cr Example: xlab = 'IVname'}

  \item{ylim}{
  \code{}(optional) The y-axis limits for the plot. 
  \cr Example: ylim = c(0, 80)}

  \item{ylab}{
  \code{}(optional) A y-axis label for the plot. 
  \cr Example: ylab = 'DVname'}

  \item{title}{
  \code{}(optional) A title for the plot. 
  \cr Example: title = 'OLS prediction of DV'}

  \item{plot_save}{
  \code{}Should a plot be saved to disk? TRUE or FALSE (the default).}

  \item{plot_save_type}{
  \code{}The output format if plot_save = TRUE. The options are 'bitmap', 'tiff', 
         'png' (the default), 'jpeg', and 'bmp'.}

  \item{cols_user}{
  \code{}A vector of colours for the levels of IV_focal_1 or IV_focal_2. 
  \cr If NULL, the default colours are selected, in order, from this vector: 
   cols_user <- c("mediumvioletred", 'black', "blue", 'cyan2', "red", 'limegreen', "yellow", 'blueviolet').
   If there are more than 7 levels of levels of IV_focal_1 or IV_focal_2, 
   then "rainbow" is used to determine the colours.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE}
}
\details{
A plot with both IV_focal_1 and IV_focal_2 predictor variables will look like an
interaction plot. But it is only a true interaction plot if the required product
term(s) was entered as a predictor when the model was created.
}
\value{A matrix with the levels of the variables that were used for the plot along 
with the predicted values, confidence intervals, and se.fit values.
}
\author{Brian P. O'Connor}
\examples{
ols_GS <- 
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               hierarchical = list( step1=c('age','quads','gluts','abdoms'), 
                                    step2=c('arms','grip')) )

PLOT_MODEL(modobject = ols_GS, 
           IV_focal_1 = 'gluts', IV_focal_1_values=NULL,
           IV_focal_2 = 'age', IV_focal_2_values=NULL, 
           IVs_nonfocal_values = NULL,
           bootstrap=TRUE, N_sims=100, CI_level=95, 
           ylim=NULL, ylab=NULL, title=NULL,
           verbose=TRUE) 
\donttest{	
ols_LW <- 
MODERATED_REGRESSION(data=data_Lorah_Wong_2018, DV='suicidal', IV='burden', MOD='belong_thwarted', 
                     IV_range='tumble',
                     MOD_levels='quantiles',
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     COVARS='depression', 
                     plot_type = 'interaction', DV_range = c(1,1.25)) 
                     
PLOT_MODEL(modobject = ols_LW, 
           IV_focal_1 = 'burden', IV_focal_1_values=NULL,
           IV_focal_2 = 'belong_thwarted', IV_focal_2_values=NULL, 
           bootstrap=TRUE, N_sims=100, CI_level=95) 
                     
logmod_Meyers <- 
  LOGISTIC_REGRESSION(data = data_Meyers_2013, DV='graduated', 
                      forced = c('sex','family_encouragement')) 

PLOT_MODEL(modobject = logmod_Meyers, 
           IV_focal_1 = 'family_encouragement', IV_focal_1_values=NULL,
           IV_focal_2=NULL, IV_focal_2_values=NULL, 
           bootstrap=FALSE, N_sims=100, CI_level=95) 
           
pois_Krem <-
  COUNT_REGRESSION(data=data_Kremelburg_2011, DV='OVRJOYED', forced=NULL, 
                   hierarchical= list(step1=c('AGE', 'female'), 
                                      step2=c('EDUC','REALRINC','DEGREE')) )

PLOT_MODEL(modobject = pois_Krem, 
           IV_focal_1 = 'AGE', 
           IV_focal_2 = 'DEGREE',
           IVs_nonfocal_values = list(EDUC = 5, female = '2'),
           bootstrap=FALSE, N_sims=100, CI_level=95) 
}
}
