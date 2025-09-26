\name{MODERATED_REGRESSION}
\alias{MODERATED_REGRESSION}
\alias{MODERATED.REGRESSION}
\title{Moderated multiple regression}
\description{Conducts moderated regression analyses for two-way interactions with
    extensive options for interaction plots, including Johnson-Neyman
    regions of significance. The output includes the
    Anova Table (Type III tests), standardized coefficients, 
    partial and semi-partial correlations, collinearity statistics,
    casewise regression diagnostics, plots of residuals and regression diagnostics, 
    and detailed information about simple slopes.
    The output includes Bayes Factors and, if requested, regression coefficients from
    Bayesian Markov Chain Monte Carlo (MCMC) analyses.}
\usage{MODERATED_REGRESSION(data, DV, IV, MOD,
                     IV_type = 'numeric', IV_range = 'tumble',
                     MOD_type='numeric', MOD_levels='quantiles', MOD_range=NULL,
                     quantiles_IV = c(.1, .9), quantiles_MOD = c(.25, .5, .75),
                     COVARS = NULL,
                     center = TRUE, 
                     CI_level = 95,
                     MCMC = FALSE,
                     Nsamples = 10000,
                     plot_type = 'residuals', plot_title=NULL, DV_range = NULL,
                     Xaxis_label = NULL, Yaxis_label=NULL, legend_label=NULL,
                     JN_type = 'Huitema', 
                     verbose = TRUE )
}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases and the columns are the variables.}

  \item{DV}{
  \code{}The name of the dependent variable. 
  \cr Example: DV = 'outcomeVar'}

  \item{IV}{
  \code{}The name of the independent variable.
  \cr Example: IV = 'varA'}
  
  \item{MOD}{
  \code{}The name of the moderator variable
  \cr Example: MOD = 'varB'}
  	    
  \item{IV_type}{
  \code{}(optional) The type of independent variable. The
   options are 'numeric' (the default) or 'factor'. 
  \cr Example: IV_type = 'factor'}
  
  \item{IV_range}{
  \code{}(optional) The independent variable range for a moderated regression plot.
  The options are:
  \itemize{
    \item{'tumble' (the default), for tumble graphs following Bodner (2016)}
    \item{'quantiles', in which case the 10th and 90th quantiles of the IV will be used
               (alternative values can be specified using the quantiles_IV argument);}
    \item{'AikenWest', in which case the IV mean - one SD, and the IV mean + one SD, will be used;}
    \item{a vector of two user-provided values (e.g., c(1, 10)); and}
    \item{NULL, in which case the minimum and maximum IV values will be used.}
  }
  Example: IV_range = 'AikenWest'}
  	    
  \item{MOD_type}{
  \code{}(optional) The type of moderator variable. The
  options are 'numeric' (the default) or 'factor'. 
  \cr Example: MOD_type = 'factor'}
  	    
  \item{MOD_levels}{
  \code{}(optional) The levels of the moderator variable to be used if MOD is continuous.
  The options are:
  \itemize{
    \item{'quantiles', in which case the .25, .5, and .75 quantiles of the MOD variable will be used
               (alternative values can be specified using the quantiles_MOD argument);}
    \item{'AikenWest', in which case the mean of MOD, the mean of MOD - one SD, and the 
               mean of MOD + one SD, will be used; and}
    \item{a vector of two user-provided values.}
    }
  Example: MOD_levels = c(1, 10)}
  	                     
  \item{MOD_range}{
  \code{}(optional) The range of the MOD values to be used in the Johnson-Neyman regions 
  of significance analyses. The options are:
  NULL (the default), in which case the minimum and maximum MOD values will be used; and
  a vector of two user-provided values. 
  \cr Example: MOD_range = c(1, 10)}
  	    
  \item{quantiles_IV}{
  \code{}(optional) The quantiles of the independent variable to be used as the IV range for 
  a moderated regression plot.
  \cr Example: quantiles_IV = c(.10, .90)}
  	    
  \item{quantiles_MOD}{
  \code{}(optional) The quantiles the moderator variable to be used as the MOD simple slope  
  values in the moderated regression analyses.
  \cr Example: quantiles_MOD = c(.25, .5, .75)}
  	    
  \item{COVARS}{
  \code{}(optional) The name(s) of possible covariates.
  \cr Example: COVARS = c('CovarA', 'CovarB', 'CovarC')}
  	    
  \item{center}{
  \code{}(optional) Logical, indicating whether the IV and MOD variables should be centered
  (default = TRUE).
  \cr Example: center = FALSE}
  	    
 \item{CI_level}{
  \code{}(optional) The confidence interval for the output, in whole numbers.
  CI_level is also used in the Johnson-Neyman regions of significance computations.
  The default is 95.}

 \item{MCMC}{
  \code{}(logical) Should Bayesian MCMC analyses be conducted? The default is FALSE.}

 \item{Nsamples}{
  \code{}(optional) The number of samples for MCMC analyses. The default is 10000.}

  \item{plot_type}{
  \code{}(optional)  The kind of plot, if any. The options are:
  \itemize{
    \item{'residuals' (the default)}
    \item{'diagnostics' (for regression diagnostics)}
    \item{'interaction' (for a traditional moderated regression interaction plot)}
    \item{'regions' (for a moderated regression Johnson-Neyman regions of significance plot), and}
    \item{'none' (for no plots).}
  }
  Example: plot_type = 'diagnostics'}

  \item{plot_title}{
  \code{}(optional) The plot title.
  \cr Example: plot_title = 'Interaction Plot'}
  
  \item{DV_range}{
  \code{}(optional) The range of Y-axis values for the plot. 
  \cr Example: DV_range = c(1,10)}
  
  \item{Xaxis_label}{
  \code{}(optional) A label for the X axis to be used in the requested plot.
  \cr Example: Xaxis_label = 'IV name'}
  
  \item{Yaxis_label}{
  \code{}(optional) A label for the Y axis to be used in the requested plot.
  \cr Example: Yaxis_label = 'DV name'}
  
  \item{legend_label}{
  \code{}(optional) A legend label for the plot.
  \cr Example: legend_label = 'MOD name'}
  
  \item{JN_type}{
  \code{}(optional) The formula to be used in computing the critical F value for the
  Johnson-Neyman regions of significance analyses. The options are 'Huitema' (the default),
  or 'Pedhazur'.
  \cr Example: JN_type = 'Pedhazur'}
  	    
  \item{verbose}{
  \code{}Should detailed results be displayed in console? The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}
}
\details{
	The Bayesian MCMC analyses can be time-consuming for larger datasets. The MCMC 
	analyses are conducted using functions, and their default settings, from the 
	BayesFactor package (Morey & Rouder, 2024).
}
\value{An object of class "MODERATED_REGRESSION". The object is a list containing the
	following possible components:
\item{modelMAINsum}{All of the summary.lm function output for the regression model 
	without interaction terms.}
\item{anova_table}{Anova Table (Type III tests).}
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
       {Bodner, T. E. (2016). Tumble graphs: Avoiding misleading end point extrapolation when  
       	graphing interactions from a moderated multiple regression analysis. 
       	\emph{Journal of Educational and Behavioral Statistics, 41,} 593-604.}
\cr\cr {Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	    multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	    Lawrence Erlbaum Associates.}
\cr\cr {Darlington, R. B., & Hayes, A. F. (2017). \emph{Regression analysis and linear models: 
	    Concepts, applications, and implementation.} Guilford Press.}
\cr\cr {Hayes, A. F. (2018a). \emph{Introduction to mediation, moderation, and conditional process 
	    analysis: A regression-based approach} (2nd ed.). Guilford Press.}
\cr\cr {Hayes, A. F., & Montoya, A. K. (2016). A tutorial on testing, visualizing, and probing 
	    an interaction involving a multicategorical variable in linear regression analysis. 
	    \emph{Communication Methods and Measures, 11,} 1-30.}
\cr\cr {Lee M. D., & Wagenmakers, E. J. (2014) \emph{Bayesian cognitive modeling: A practical 
        course.} Cambridge University Press.}
\cr\cr {Morey, R. & Rouder, J. (2024). \emph{BayesFactor: Computation of Bayes Factors for 
        Common Designs.} R package version 0.9.12-4.7, 
        https://github.com/richarddmorey/bayesfactor.}
\cr\cr {O'Connor, B. P. (1998). All-in-one programs for exploring interactions in moderated 
	multiple regression. \emph{Educational and Psychological Measurement, 58,} 833-837.}
\cr\cr {Pedhazur, E. J. (1997). \emph{Multiple regression in behavioral research: Explanation 
   and prediction.} (3rd ed.). Wadsworth Thomson Learning.}
}
\author{Brian P. O'Connor}
\examples{
# moderated regression	-- with IV_range = 'AikenWest'
MODERATED_REGRESSION(data=data_Lorah_Wong_2018, DV='suicidal', IV='burden',  MOD='belong_thwarted', 
                     IV_range='AikenWest',
                     MOD_levels='quantiles',
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = TRUE, COVARS='depression', 
                     plot_type = 'interaction', plot_title=NULL, DV_range = c(1,1.25))
\donttest{
# moderated regression	-- with  IV_range = 'tumble'
MODERATED_REGRESSION(data=data_Lorah_Wong_2018, DV='suicidal', IV='burden', MOD='belong_thwarted', 
                     IV_range='tumble',
                     MOD_levels='quantiles',
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = TRUE, COVARS='depression', 
                     plot_type = 'interaction', plot_title=NULL, DV_range = c(1,1.25)) 

# moderated regression	-- with numeric values for IV_range & MOD_levels='AikenWest'       
MODERATED_REGRESSION(data=data_OConnor_Dvorak_2001, DV='Aggressive_Behavior', 
                     IV='Maternal_Harshness', MOD='Resiliency', 
                     IV_range=c(1,7.7), 
                     MOD_levels='AikenWest', MOD_range=NULL,
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = FALSE, 
                     plot_type = 'interaction', 
                     DV_range = c(1,6), 
                     Xaxis_label='Maternal Harshness', 
                     Yaxis_label='Adolescent Aggressive Behavior', 
                     legend_label='Resiliency')
}
}
