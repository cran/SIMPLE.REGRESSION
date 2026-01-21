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
\usage{OLS_REGRESSION(data, DV, forced=NULL, hierarchical=NULL, formula=NULL, 
               CI_level = 95,
               MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                   thin = 1, burnin = 1000, 
                                   HDI_plot_est_type = 'standardized'),
               plot_type = 'residuals', 
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

  \item{formula}{
  \code{}(optional) Text for an R formula. Useful for testing for interactions.
  \cr Example: formula = "Aggressive_Behavior ~ Maternal_Harshness * Resiliency"")}

 \item{CI_level}{
  \code{}(optional) The confidence interval for the output, in whole numbers.
  The default is 95.}

 \item{MCMC_options}{
  \code{}(optional) A list specifying the following options for Bayesian MCMC
  analyses: (1) "MCMC", Should MCMC analyses be conducted? The options are TRUE
  or FALSE; (2) "Nsamples", for the number of iterations or samples from the
  posterior distribution; (3) "thin", for the chain-thinning interval;
  (4) "burnin", for the burnin period, i.e., the number of initial samples
  that should be dropped from the chains; and (5) "HDI_plot_est_type", for
  the kind of regression estimates that will appear in any requested HDI 
  plots. The options are "standardized" or "raw".
  \cr Example: MCMC_options = list(MCMC = TRUE, Nsamples = 10000, thin = 1, 
                                   burnin = 1000, 
                                   HDI_plot_est_type = 'standardized')}

  \item{plot_type}{
  \code{}(optional)  The kind of plots, if any. The options are:
  \itemize{
    \item{'residuals' (the default)}
    \item{'diagnostics' (for regression diagnostics),}
    \item{'Bayes_HDI' (for MCMC posterior distributions), or}
    \item{'none' (for no plots).}
  }
  Example: plot_type = 'diagnostics'}

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
	
	The function assigns contrasts (dummy codes) to factor variables that do not
	already have contrasts. The baseline group for the dummy codes is determined 
	by the alphabetic/numeric order of the factor levels. If the terms "control"
	or "Control" or "baseline" or "Baseline" appear in the names of a factor
	level, then that factor level is used as the dummy codes baseline.
	
	The MCMC analyses are conducted using functions, and their default settings,  
	from the rstanarm	package (Goodrich, Gabry, Ali, & Brilleman, 2024).  
	MCMC analyses can be time-consuming for larger datasets.
	
	The Bayes Factor analyses are conducted using functions, and their default  
	settings, from the BayesFactor package (Morey & Rouder, 2024). 
	
  The Bayes factor values for the predictor variables are based solely on the
  predictor t & df values, using the ttest.tstat function from the BayesFactor
  package.

  Good sources for interpreting residuals and diagnostics plots:
	\itemize{
  \item{\href{https://library.virginia.edu/data/articles/diagnostic-plots}{library.virginia.edu}}
  \item{\href{ https://www.sthda.com/english/articles/index.php?url=/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/}{sthda.com}}
  }
}
\value{An object of class "OLS_REGRESSION". The object is a list containing the
	following possible components:
\item{model}{All of the lm function output for the regression model 
	without interaction terms.}
\item{modelsum}{All of the summary.lm function output for the regression model 
	without interaction terms.}
\item{anova_table}{Anova Table (Type III tests).}
\item{mainRcoefs}{Predictor coefficients for the model without interaction terms.}
\item{modeldata}{All of the predictor and outcome raw data that were used in the model,
	along with regression diagnostic statistics for each case.}
\item{collin_diags}{Collinearity diagnostic coefficients for models without interaction terms.}
\item{chain_dat}{The MCMC chains.}
\item{Bayes_HDIs}{The Bayesian HDIs.}
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
\cr\cr {Goodrich, B., Gabry, J., Ali, I., & Brilleman, S. (2024). \emph{rstanarm:  
        Bayesian applied regression modeling via Stan.} R package version 2.32.1, 
        https://mc-stan.org/rstanarm/.}
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
# for Green_Salkind_2014, p. 263: forced (simultaneous) entry
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               forced = c('quads','gluts','abdoms','arms','grip'))
\donttest{
# for Green_Salkind_2014, p. 265: hierarchical entry with Bayesian MCMC analyses & HDI plots
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               hierarchical = list(step1=c('quads','gluts','abdoms'), 
                                   step2=c('arms','grip')),
               MCMC_options = list(MCMC = TRUE, Nsamples = 10000, 
                                   thin = 1, burnin = 1000, 
                                   HDI_plot_est_type = 'raw'),
               plot_type = 'Bayes_HDI')
               
# for O'Connor & Dvorak, 2001, p. 17; 2-way interaction specified via formula      
OLS_REGRESSION(data=data_OConnor_Dvorak_2001, 
               formula = 'Aggressive_Behavior ~ Maternal_Harshness * Resiliency')
}
}
