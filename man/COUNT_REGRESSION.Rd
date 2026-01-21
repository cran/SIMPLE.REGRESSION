\name{COUNT_REGRESSION}
\alias{COUNT_REGRESSION}
\title{Count data regression}
\description{Provides SPSS- and SAS-like output for count data
   regression, including Poisson, quasi-Poisson, negative binomial,
   zero-inflated poisson, zero-inflated negative binomial,
   hurdle poisson, and hurdle negative binomial models.  
   The output includes model summaries, classification tables, omnibus tests of 
   the model coefficients, overdispersion tests, model effect sizes, the 
   model coefficients, correlation matrix for the model coefficients, collinearity 
   statistics, and casewise regression diagnostics.}
\usage{COUNT_REGRESSION(data, DV, forced = NULL, hierarchical = NULL, formula=NULL,
                 model_type = 'poisson',
                 offset = NULL,
                 CI_level = 95,
                 MCMC_options = list(MCMC = FALSE, Nsamples = 10000, 
                                     thin = 1, burnin = 1000, 
                                     HDI_plot_est_type = 'raw'),
                 plot_type = 'residuals',
                 GoF_model_types = TRUE,
                 verbose = TRUE )
}
\arguments{
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

  \item{model_type}{
  \code{}(optional) The name of the error distribution to be 
  used in the model. The options are:
  \itemize{
   \item{"poisson" (the default),}
   \item{"quasipoisson",}
   \item{"negbin", for negative binomial,}
   \item{"zinfl_poisson", for zero-inflated poisson,}
   \item{"zinfl_negbin", for zero-inflated negative binomial, or}
   \item{"hurdle_poisson", for hurdle poisson, or}
   \item{"hurdle_negbin", for hurdle negative binomial.}}

  Example: model_type = 'quasipoisson'}

  \item{offset}{
  \code{}(optional) The name of the offset variable, if there is one. This variable
  should be in the desired metric (e.g., log). No transformation of an
  offset variable is performed internally.
  \cr Example: offset = 'Varname'}

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

 \item{GoF_model_types}{
  \code{}(optional) Should fit coefficients be computed for multiple model 
  types (Poisdon, quasi-Poisson, negative binomial, zero-inflated Poisson, 
  zero-inflated negative binomial, and hurdle)? The default is TRUE.}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}
}
\details{
	This function uses the glm function from the stats package, the 
	negative.binomial function from the MASS package, and the zeroinfl
	and hurdle functions from the pscl package (Zeileis, Kleiber, & Jackman, 2008).
  It supplements the output from these packages with additional statistics and 
	in formats that resemble SPSS and SAS output. The predictor variables can 
	be numeric or factors.
	
	The function assigns contrasts (dummy codes) to factor variables that do not
	already have contrasts. The baseline group for the dummy codes is determined 
	by the alphabetic/numeric order of the factor levels. If the terms "control"
	or "Control" or "baseline" or "Baseline" appear in the names of a factor
	level, then that factor level is used as the dummy codes baseline.
	
  The following descriptions of zero-inflated and hurdle models were
  provided by Atkins and Baldwin (2013), by Friendly and Meyer (2016), and at
  https://stats.oarc.ucla.edu/r/dae/zinb/:

  Zero-inflated and hurdle models are used when there is an overabundance of 
  zero counts (excessive, or over-dispersed zero counts). Both have two
  submodels, one related to the zeroes and a second related to the
  counts. The key difference between hurdle and zero-inflated models is
  how they handle zeroes: Hurdle models cleanly divide the models, with
  all zeroes accounted for in the logistic regression, whereas
  zero-inflated models treat the observed zeroes as a mixture from two latent 
  classes that produce zeroes.

  \strong{Zero-inflated models} assume that the observed counts arise from a
  mixture of two latent classes of observations: some structural zeros for
  whom the DV will always be zero, and a second class for whom the observed
  count may be zero or above zero. The excess zeros are assumed to have
  been generated by a separate process from the count values and it is
  assumed that the excess zeros can be modeled independently.

  For example, imagine that wildlife biologists want to model how many
  fish are being caught by visitors to a park. Some visitors do not fish
  (structural zeros), but there is no data on whether a person fished or
  not. Some visitors who did fish did not catch any fish so there are
  excess zeros in the data because of the people that did not fish. The
  variables that predict whether or not visitors fished may or may not be
  the same variables that predict how many fish visitors caught. Separate
  models for the zeroes and for the counts can be examined. Zero-inflated
  models assume that zero values are due to two dierent processes, e.g.,
  that a visitor has gone fishing vs. not gone fishing. If not gone
  fishing, the only outcome possible is zero. If gone fishing, it is then
  a count process. The two parts of the a zero-inflated model are a binary
  (logistic) model and a count model (Poisson or negative binomial). The
  expected counts are expressed as a combination of the two processes.

	For the zero (logistic) portion of zero-inflated models, the predicted outcomes
	are the zero values (excess zeros) for the DV. A positive coefficient (B) for 
	a predictor thus means that as values on a predictor increase, the probability 
	of observing a zero value for the DV increases.

  \strong{Hurdle models} also deal with an excess of zero DV values, but without
  assuming that zero values arise from a mixture of two latent classes of
  observations. Imagine that it is (somehow) known that every visitor to a
  park did in fact fish. There could be an excess of zeroes because many
  of the visitors did not know how to fish. A separate logistic regression
  submodel is used to distinguish zero counts from the larger counts. The
  submodel for the positive counts is a truncated Poisson or
  negative-binomial model, excluding the zero counts. In other words,
  there is one process and submodel accounting for the zero counts and a
  separate process accounting for the positive counts, once the zero hurdle
  has been crossed. In zero-inflation models, the first process
  generates only extra zeros beyond those of the regular Poisson
  distribution. For hurdle models, the first process generates all of the
  zeros. In hurdle models, the zero values are considered fully observed,
  rather than latent.

	For the zero (logistic) portion of hurdle models, the predicted outcomes
	are for going from zero to greater than zero values for the DV. A positive 
	coefficient (B) for a predictor thus means that as values on a predictor increase,
  the probability of crossing the hurdle (obtaining a value higher than zero) 
  for the DV increases.
	
	\strong{Predicted values}, for selected levels of the predictor variables,
	can be produced and plotted using the PLOT_MODEL funtion in this package.

	\strong{The Bayesian MCMC analyses} can be time-consuming for larger datasets. The MCMC 
	analyses are conducted using functions, and their default settings, from the 
	rstanarm package (Goodrich, Gabry, Ali, & Brilleman, 2024). 
	Family = 'quasipoisson' analyses are currently not possible for the MCMC 
	analyses. model_type = 'poisson' is therefore used instead.
	
	The Bayesian MCMC analyses are currently not available for zero-inflated 
	poisson and zero-inflated negative binomial models.
	
	The MCMC results can be verified using the model checking functions in the
	rstanarm package (e.g., Muth, Oravecz, & Gabry, 2018).

  Good sources for interpreting count data regression residuals and diagnostics plots:
	\itemize{
  \item{\href{https://rpubs.com/benhorvath/glm_diagnostics}{rpubs.com/benhorvath}}
  \item{\href{https://library.virginia.edu/data/articles/understanding-deviance-residuals}{library.virginia.edu}}
  }
}
\value{An object of class "COUNT_REGRESSION". The object is a list containing the
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
       {Atkins, D. C., Baldwin, S. A., Zheng, C., Gallop, R. J., & Neighbors, C. (2013). 
       A tutorial on count regression and zero-altered count models for 
       longitudinal substance use data. \emph{Psychology of Addictive Behaviors, 27(1),}
       166177. https://doi.org/10.1037/a0029508}
       \cr\cr {Atkins, D. C., & Gallop, R. J. (2007). Rethinking how family researchers 
        model infrequent outcomes: A tutorial on count regression and zero-inflated 
        models. \emph{Journal of Family Psychology, 21(4),} 726-735.}
\cr\cr {Beaujean, A. A., & Grant, M. B. (2019). Tutorial on using regression 
        models with count outcomes using R. \emph{Practical Assessment, 
        Research, and Evaluation: Vol. 21, Article 2.}}
\cr\cr {Coxe, S., West, S.G., & Aiken, L.S. (2009). The analysis of count data: 
        A gentle introduction to Poisson regression and its alternatives. 
       	\emph{Journal of Personality Assessment, 91,} 121-136.}
\cr\cr {Dunn, P. K., & Smyth, G. K. (2018). \emph{Generalized linear models 
       with examples in R.} Springer.}
\cr\cr {Friendly, M., & Meyer, D. (2016). \emph{Discrete Data Analysis with R:  
       Visualization and Modeling Techniques for Categorical and Count Data.} 
       Chapman and Hall/CRC.}
\cr\cr {Hardin, J. W., & Hilbe, J. M. (2007). \emph{Generalized linear models 
        and extensions.} Stata Press.}
\cr\cr {Muth, C., Oravecz, Z., & Gabry, J. (2018). User-friendly Bayesian regression 
        modeling: A tutorial with rstanarm and shinystan. \emph{The Quantitative Methods 
        for Psychology, 14(2),} 99119. 
        \cr https://doi.org/10.20982/tqmp.14.2.p099}
\cr\cr {Orme, J. G., & Combs-Orme, T. (2009). \emph{Multiple regression with discrete 
        dependent variables.} Oxford University Press.}
\cr\cr {Rindskopf, D. (2023). Generalized linear models. In H. Cooper, M. N. 
        Coutanche, L. M. McMullen, A. T. Panter, D. Rindskopf, & K. J. Sher (Eds.), 
	     \emph{APA handbook of research methods in psychology: Data analysis and 
	     research publication, } (2nd ed., pp. 201-218). American Psychological Association.}
\cr\cr {Zeileis, A., Kleiber, C., & Jackman, S. (2008). Regression Models for Count Data in R. 
        \emph{Journal of Statistical Software, 27(8).} https://www.jstatsoft.org/v27/i08/.}
}
\author{Brian P. O'Connor}
\examples{
# for Kremelburg, 2011, p.262 & p. 282: poisson regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='OVRJOYED', 
                 forced=c('AGE','EDUC','REALRINC','female'))

\donttest{
# for Kremelburg, 2011, p. 266-267 & p. 284: 
# negative binomial regression & with Bayesian MCMC analyses & HDI plots
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='HURTATWK', 
                 forced=c('AGE','EDUC','REALRINC','female'),
                 model_type = 'negbin',
                 MCMC_options = list(MCMC = TRUE, Nsamples = 10000, 
                                     thin = 1, burnin = 1000, 
                                     HDI_plot_est_type = 'raw'),
                 plot_type = 'Bayes_HDI')

# for Orme, 2009, p. 160: with an offset variable
COUNT_REGRESSION(data=data_Orme_2009_5, DV='NumberAdopted', forced=c('Married'), 
                 offset='lnYearsFostered')

# zero-inflated poisson regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='HURTATWK', 
                 forced=c('AGE','EDUC','REALRINC','female'),
                 model_type = 'zinfl_poisson',
                 plot_type = 'diagnostics')

# hurdle negative binomial regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='HURTATWK', 
                 forced=c('AGE','EDUC','REALRINC','female'),
                 model_type = 'hurdle_negbin',
                 plot_type = 'diagnostics')
}
}
