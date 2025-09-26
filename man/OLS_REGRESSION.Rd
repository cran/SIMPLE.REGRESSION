\name{OLS_REGRESSION}
\alias{OLS_REGRESSION}
\alias{SIMPLE.REGRESSION}
\title{Ordinary least squares regression}
\description{Provides SPSS- and SAS-like output for ordinary least squares simultaneous 
	entry regression and hierarchical entry regression. The output includes the
  Anova Table (Type III tests), standardized coefficients, 
  partial and semi-partial correlations, collinearity statistics,
  casewise regression diagnostics, plots of residuals and regression diagnostics.
  The output includes Bayes Factors and, if requested, regression coefficients from
  Bayesian Markov Chain Monte Carlo (MCMC) analyses.}
\usage{OLS_REGRESSION(data, DV, forced=NULL, hierarchical=NULL, 
               COVARS=NULL,
               plot_type = 'residuals', 
               CI_level = 95,
               MCMC = FALSE,
               Nsamples = 10000,
               verbose=TRUE, \dots)
}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases and the columns are the variables.}

  \item{DV}{
  \code{}The name of the dependent variable. 
  \cr Example: DV = 'outcomeVar'}

  \item{forced}{
  \code{}(optional) A vector of the names of the predictor variables for a forced/simultaneous  
  entry regression. The variables can be numeric or factors.
  \cr Example: forced = c('VarA', 'VarB', 'VarC')}

  \item{hierarchical}{
  \code{}(optional) A list with the names of the predictor variables for each step of 
  a hierarchical regression. The variables can be numeric or factors.
  \cr Example: hierarchical = list(step1=c('VarA', 'VarB'), step2=c('VarC', 'VarD'))}

  \item{COVARS}{
  \code{}(optional) The name(s) of possible covariates variable for a moderated regression
  analysis.
  \cr Example: COVARS = c('CovarA', 'CovarB', 'CovarC')}
  	    
  \item{plot_type}{
  \code{}(optional)  The kind of plots, if any. The options are:
  \itemize{
    \item{'residuals' (the default)}
    \item{'diagnostics' (for regression diagnostics), or}
    \item{'none' (for no plots).}
  }
  Example: plot_type = 'diagnostics'}

 \item{CI_level}{
  \code{}(optional) The confidence interval for the output, in whole numbers.
  The default is 95.}

 \item{MCMC}{
  \code{}(logical) Should Bayesian MCMC analyses be conducted? The default is FALSE.}

 \item{Nsamples}{
  \code{}(optional) The number of samples for MCMC analyses. The default is 10000.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}

  \item{\dots}{(dots, for internal purposes only at this time.)}
}
\details{
	This function uses the lm function from the stats package,
	supplements the output with additional statistics, and it formats the output
	so that it resembles SPSS and SAS regression output. The predictor
	variables can be numeric or factors.
	
	The Bayesian MCMC analyses can be time-consuming for larger datasets. The MCMC 
	analyses are conducted using functions, and their default settings, from the 
	BayesFactor	package (Morey & Rouder, 2024). 

  Good sources for interpreting residuals and diagnostics plots:
	\itemize{
  \item{\href{https://library.virginia.edu/data/articles/diagnostic-plots}{library.virginia.edu}}
  \item{\href{ https://www.sthda.com/english/articles/index.php?url=/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/}{sthda.com}}
  }
}
\value{An object of class "OLS_REGRESSION". The object is a list containing the
	following possible components:
\item{modelMAIN}{All of the lm function output for the regression model 
	without interaction terms.}
\item{modelMAINsum}{All of the summary.lm function output for the regression model 
	without interaction terms.}
\item{anova_table}{Anova Table (Type III tests).}
\item{mainRcoefs}{Predictor coefficients for the model without interaction terms.}
\item{modeldata}{All of the predictor and outcome raw data that were used in the model,
	along with regression diagnostic statistics for each case.}
\item{collin_diags}{Collinearity diagnostic coefficients for models without interaction terms.}
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
# forced (simultaneous) entry
head(data_Green_Salkind_2014)
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               forced = c('quads','gluts','abdoms','arms','grip'))
\donttest{
# hierarchical entry
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               hierarchical = list( step1=c('quads','gluts','abdoms'), 
                                    step2=c('arms','grip')) )
}
}
