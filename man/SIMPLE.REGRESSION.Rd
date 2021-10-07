 \name{SIMPLE.REGRESSION}
\alias{SIMPLE.REGRESSION}
\title{Multiple regression and moderated multiple regression}
\description{Provides SPSS- and SAS-like output for least squares simultaneous 
	entry regression, hierarchical entry regression, and
    moderated regression, as well as interaction plots and Johnson-Neyman
    regions of significance for interactions. The output includes standardized
    coefficients, partial and semi-partial correlations, collinearity diagnostics,
    plots of residuals, and detailed information about simple slopes for interactions.}
\usage{SIMPLE.REGRESSION(data, DV, forced=NULL, hierarchical=NULL, 
        IV=NULL, IV_type = 'numeric', IV_range='tumble',
        MOD=NULL, MOD_type = 'numeric', MOD_levels='quantiles', MOD_range=NULL,
        quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
        CENTER = TRUE, COVARS=NULL,
        PLOT_type = 'regions', PLOT_title=NULL,
        Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
        JN_type = 'Huitema', verbose=TRUE )
}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases and the columns are the variables.}

  \item{DV}{
  \code{}The name of the dependent variable, e.g., DV = 'outcomeVar'.}

  \item{forced}{
  \code{}(optional) A vector of the names of the predictor variables for a forced/simultaneous  
  entry regression, e.g., forced = c('VarA', 'VarB', 'VarC').
  \cr The variables can be numeric or factors.}

  \item{hierarchical}{
  \code{}(optional) A list with the names of the predictor variables for each step of a hierarchical regression,
  e.g., hierarchical = \cr  list(step1=c('VarA', 'VarB'), step2=c('VarC', 'VarD'),  step3=c('VarE', 'VarF')).
  \cr The variables can be numeric or factors.}

  \item{IV}{
  \code{}(optional) The name of the independent variable for a moderated regression. Not required 
  for forced or hierarchical regression.}
  	    
  \item{IV_type}{
  \code{}(optional) The type of independent variable for a moderated regression. The
  options are 'numeric' (the default) or 'factor'. Not required 
  for forced or hierarchical regression.}

  \item{IV_range}{
  \code{}(optional) The independent variable range for a moderated regression plot.
  The options are:
  'tumble' (the default), for tumble graphs following Bodner (2016);
  'quantiles', in which case the 10th and 90th quantiles of the IV will be used
               (alternative values can be specified using the quantiles_IV argument);
  NULL, in which case the minimum and maximum IV values will be used;
  'AikenWest', in which case the IV mean - one SD, and the IV mean + one SD, will be used; and
  a vector of two user-provided values (e.g., c(1, 10)). 
  }
  	    
  \item{MOD}{
  \code{}(optional) The name of the moderator variable for a moderated regression. Not required 
  for a regular (non moderated) multiple regression.}
  	    
  \item{MOD_type}{
  \code{}(optional) The type of moderator variable for a moderated regression. The
  options are 'numeric' (the default) or 'factor'. Not required 
  for forced or hierarchical regression.}
  	    
  \item{MOD_levels}{
  \code{}(optional) The levels of the moderator variable to be used in a moderated regression,
  if MOD is continuous. Not required for a regular (non moderated) multiple regression.
  The options are:
  'quantiles', in which case the .25, .5, and .75 quantiles of the MOD variable will be used
               (alternative values can be specified using the quantiles_MOD argument);
  'AikenWest', in which case the mean of MOD, the mean of MOD - one SD, and the 
               mean of MOD + one SD, will be used; and
  a vector of two user-provided values (e.g., c(1, 10)). 
  }
  	                     
  \item{MOD_range}{
  \code{}(optional)  The range of the MOD values to be used in the Johnson-Neyman regions 
  of significance analyses. The options are:
  NULL (the default), in which case the minimum and maximum MOD values will be used; and
  a vector of two user-provided values (e.g., c(1, 10)). 
  }
  	    
  \item{quantiles_IV}{
  \code{}(optional)  The quantiles the independent variable to be used as the IV range for 
  a moderated regression plot.}
  	    
  \item{quantiles_MOD}{
  \code{}(optional)  The quantiles the moderator variable to be used as the MOD simple slope  
  values in the moderated regression analyses.}
  	    
  \item{CENTER}{
  \code{}(optional)  Logical indicating whether the IV and MOD variables should be centered
  in a moderated regression analysis (default = TRUE).}
  	    
  \item{COVARS}{
  \code{}(optional)  The name(s) of possible covariates variable for a moderated regression
  analysis, e.g., \cr COVARS = c('CovarA', 'CovarB', 'CovarC').}
  	    
  \item{PLOT_type}{
  \code{}(optional)  The kind of plot for a moderated regression. The options are 'interaction'
  for a traditional interaction plot, 
  and 'regions', for a Johnson-Neyman regions of significance plot.}

  \item{PLOT_title}{
  \code{}(optional)  The plot title for a moderated regression.}
  
  \item{Xaxis_label}{
  \code{}(optional)  A label for the X axis to be used in the requested plot.}
  
  \item{Yaxis_label}{
  \code{}(optional)  A label for the Y axis to be used in the requested plot.}
  
  \item{LEGEND_label}{
  \code{}(optional)  The legend label for a moderated regression.}
  
  \item{JN_type}{
  \code{}(optional)  The formula to be used in computing the critical F value for the
  Johnson-Neyman regions of significance analyses. The options are 'Huitema' (the default),
  or 'Pedhazur'.}
  	    
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}
}
\details{
	This function relies heavily on the lm function from the stats package.
	It supplements the lm function output with additional statistics and it formats the output
	so that it resembles SPSS and SAS regression output. The predictor
	variables can be numeric or factors. Only least squares regressions are performed.
}
\value{An object of class "SIMPLE.REGRESSION". The object is a list containing the
	following possible components:
\item{modelMAIN}{All of the lm function output for the regression model 
	without interaction terms.}
\item{modelMAINsum}{All of the summary.lm function output for the regression model 
	without interaction terms.}
\item{mainRcoefs}{Predictor coefficients for the model without interaction terms.}
\item{modeldata}{All of the predictor and outcome raw data that were used in the model,
	along with regression diagnostic statistics for each case.}
\item{collin_diags}{Collinearity diagnostic coefficients for models without interaction terms.}
\item{modelXNsum}{Regression model statistics with interaction terms.}
\item{RsqchXn}{Rsquared change for the interaction.}
\item{fsquaredXN}{fsquared change for the interaction.}
\item{xnRcoefs}{Predictor coefficients for the model with interaction terms.}
\item{simslop}{The simple slopes.}
\item{simslopZ}{The standardized simple slopes.}
\item{plotdon}{The plot data for a moderated regression.}
\item{JN.data}{The Johnson-Neyman results for a moderated regression.}
\item{ros}{The Johnson-Neyman regions of significance for a moderated regression.}
}
\references{
       {Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	    multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	    Lawrence Erlbaum Associates Publishers.}
\cr\cr {Darlington, R. B., & Hayes, A. F. (2017). \emph{Regression analysis and linear models: 
	    Concepts, applications, and implementation.} New York: The Guilford Press.}
\cr\cr {Hayes, A. F. (2018a). \emph{Introduction to mediation, moderation, and conditional process 
	    analysis: A regression-based approach} (2nd ed.). New York, NY: Guilford Press.}
\cr\cr {Hayes, A. F., & Montoya, A. K. (2016). A tutorial on testing, visualizing, and probing 
	    an interaction involving a multicategorical variable in linear regression analysis. 
	    \emph{Communication Methods and Measures, 11,} 1-30.}
\cr\cr {Pedhazur, E. J. (1997). \emph{Multiple regression in behavioral research: Explanation 
   and prediction.} (3rd ed.). Fort Worth, Texas: Wadsworth Thomson Learning}
}
\author{Brian P. O'Connor }
\examples{
# forced (simultaneous) entry
head(data_Green_Salkind_2008)
SIMPLE.REGRESSION(data=data_Green_Salkind_2008, DV='injury', 
                  forced = c('quads','gluts','abdoms','arms','grip'))

\donttest{	
# hierarchical entry
SIMPLE.REGRESSION(data=data_Green_Salkind_2008, DV='injury', 
                  hierarchical = list( step1=c('quads','gluts','abdoms'), step2=c('arms','grip')) )
		

# moderated regression	
head(data_Lorah_Wong_2018)
SIMPLE.REGRESSION(data=data_Lorah_Wong_2018, DV='sis',
       IV='burden', IV_type = 'numeric', IV_range='tumble',
       MOD='belong', MOD_type = 'numeric', MOD_levels='quantiles', MOD_range=NULL,
       quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
       CENTER = TRUE, COVARS='dep', 
       PLOT_type = 'regions', PLOT_title=NULL, 
       Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
       JN_type = 'Huitema', verbose=TRUE ) 
}
}
