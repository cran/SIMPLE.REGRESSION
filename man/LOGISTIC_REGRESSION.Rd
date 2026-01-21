\name{LOGISTIC_REGRESSION}
\alias{LOGISTIC_REGRESSION}
\title{Logistic regression}
\description{Logistic regression analyses with SPSS- and SAS-like output. The output includes 
   model summaries, classification tables, omnibus tests of model coefficients, 
   the model coefficients, likelihood ratio tests for the predictors, overdispersion 
   tests, model effect sizes, the correlation matrix for the model coefficients,
   collinearity statistics, and casewise regression diagnostics.}
\usage{LOGISTIC_REGRESSION(data, DV, forced = NULL, hierarchical = NULL, formula=NULL,
                    ref_category = NULL,
                    family = 'binomial',
                    CI_level = 95,
                    MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                        thin = 1, burnin = 1000, 
                                        HDI_plot_est_type = 'standardized'),
                    plot_type = 'residuals',
                    verbose = TRUE)
}\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases and the columns are the variables.}

  \item{DV}{
  \code{}The name of the dependent variable. 
  \cr Example: DV = 'outcomeVar'.}

  \item{forced}{
  \code{}(optional) A vector of the names of the predictor variables for a forced/simultaneous  
  entry regression. The variables can be numeric or factors. 
  \cr Example: forced = c('VarA', 'VarB', 'VarC')}

  \item{hierarchical}{
  \code{}(optional) A list with the names of the predictor variables for each step of a 
  hierarchical regression. The variables can be numeric or factors.
  \cr Example: hierarchical = list(step1=c('VarA', 'VarB'), step2=c('VarC', 'VarD'))}

  \item{formula}{
  \code{}(optional) Text for an R formula. Useful for testing for interactions.
  \cr Example: formula = "Aggressive_Behavior ~ Maternal_Harshness * Resiliency"")}

  \item{ref_category}{
  \code{}(optional) The reference category for DV. 
  \cr Example: ref_category = 'alive'}

  \item{family}{
  \code{}(optional) The name of the error distribution to be used in the model. The options are:
  \itemize{
   \item{"binomial" (the default), or}
   \item{"quasibinomial", which should be used when there is overdispersion.}
  }
  Example: family = 'quasibinomial'}
  
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
  \code{}(optional) The kind of plots, if any. The options are:
  \itemize{
   \item{'residuals' (the default),}
   \item{'diagnostics', for regression diagnostics,}
    \item{'Bayes_HDI' (for MCMC posterior distributions), and}
   \item{'none', for no plots.}
  }
  Example: plot_type = 'diagnostics'}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}
}
\details{
	This function uses the glm function from the stats package and
	supplements the output with additional statistics and in formats that resembles
	SPSS and SAS output. The predictor variables can be numeric or factors.
	
  The function assigns contrasts (dummy codes) to factor variables that do not
	already have contrasts. The baseline group for the dummy codes is determined 
	by the alphabetic/numeric order of the factor levels. If the terms "control"
	or "Control" or "baseline" or "Baseline" appear in the names of a factor
	level, then that factor level is used as the dummy codes baseline.

	Predicted values for this model, for selected levels of the predictor variables,
	can be produced and plotted using the PLOT_MODEL funtion in this package.

	The Bayesian MCMC analyses can be time-consuming for larger datasets. The MCMC 
	analyses are conducted using functions, and their default settings, from the 
	rstanarm package (Goodrich, Gabry, Ali, & Brilleman, 2024). 
	The MCMC results can be verified using the model checking functions in the
	rstanarm package (e.g., Muth, Oravecz, & Gabry, 201).

  Good sources for interpreting logistic regression residuals and diagnostics plots:
	\itemize{
  \item{\href{https://rpubs.com/benhorvath/glm_diagnostics}{rpubs.com/benhorvath}}
  \item{\href{https://library.virginia.edu/data/articles/understanding-deviance-residuals}{library.virginia.edu}}
  }
}
\value{An object of class "LOGISTIC_REGRESSION". The object is a list containing the
	following possible components:
\item{model}{All of the glm function output for the regression model.}
\item{modelsum}{All of the summary.glm function output for the regression model.}
\item{modeldata}{All of the predictor and outcome raw data that were used in the model,
	along with regression diagnostic statistics for each case.}
\item{collin_diags}{Collinearity diagnostic coefficients for models without interaction terms.}
\item{chain_dat}{The MCMC chains.}
\item{Bayes_HDIs}{The Bayesian HDIs.}
}
\references{
       {Dunn, P. K., & Smyth, G. K. (2018). \emph{Generalized linear models 
       with examples in R.} Springer.}
\cr\cr {Field, A., Miles, J., & Field, Z. (2012). 
       	\emph{Discovering statistics using R.} Los Angeles, CA: Sage.}
\cr\cr {Goodrich, B., Gabry, J., Ali, I., & Brilleman, S. (2024). \emph{rstanarm:  
        Bayesian applied regression modeling via Stan.} R package version 2.32.1, 
        https://mc-stan.org/rstanarm/.}
\cr\cr {Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2014). 
        \emph{Multivariate data analysis,} (8th ed.).
	      Lawrence Erlbaum Associates.}
\cr\cr {Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013) 
        \emph{Applied logistic regression.} (3rd ed.). John Wiley & Sons.}
\cr\cr {Muth, C., Oravecz, Z., & Gabry, J. (2018). User-friendly Bayesian regression 
        modeling: A tutorial with rstanarm and shinystan. \emph{The Quantitative Methods 
        for Psychology, 14(2),} 99119. 
        \cr https://doi.org/10.20982/tqmp.14.2.p099}
\cr\cr {Orme, J. G., & Combs-Orme, T. (2009). \emph{Multiple regression with discrete 
        dependent variables.} Oxford University Press.}
\cr\cr {Pituch, K. A., & Stevens, J. P. (2016). 
	      \emph{Applied multivariate statistics for the social sciences: Analyses with 
      	SAS and IBM's SPSS,} (6th ed.). Routledge.}
\cr\cr {Rindskopf, D. (2023). Generalized linear models. In H. Cooper, M. N. 
        Coutanche, L. M. McMullen, A. T. Panter, D. Rindskopf, & K. J. Sher (Eds.), 
	      \emph{APA handbook of research methods in psychology: Data analysis and 
	      research publication, } (2nd ed., pp. 201-218). American Psychological Association.}
}
\author{Brian P. O'Connor}
\examples{
# Meyers, 2013, p. 263: forced (simultaneous) entry
LOGISTIC_REGRESSION(data = data_Meyers_2013, DV='graduated', 
                    forced=c('sex','family_encouragement'),
                    plot_type = 'diagnostics')
\donttest{	
# for Kremelburg, 2011, p. 244: hierarchical entry, with Bayesian MCMC analyses & HDI plots
LOGISTIC_REGRESSION(data = data_Kremelburg_2011, DV='OCCTRAIN',
                    hierarchical=list( step1=c('AGE', 'female'), 
                                       step2=c('EDUC','REALRINC')),
                    MCMC_options = list(MCMC = TRUE, Nsamples = 10000, 
                                        thin = 1, burnin = 1000, 
                                        HDI_plot_est_type = 'raw'),
                    plot_type = 'Bayes_HDI')
}
}
