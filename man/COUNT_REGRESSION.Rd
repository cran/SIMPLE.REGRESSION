\name{COUNT_REGRESSION}
\alias{COUNT_REGRESSION}
\title{Count data regression}
\description{Provides SPSS- and SAS-like output for count data
   regression, including Poisson, quasi-Poisson, and negative binomial models.  
   The output includes model summaries, classification tables, omnibus tests of 
   the model coefficients, overdispersion tests, model effect sizes, the 
   model coefficients, correlation matrix for the model coefficients, collinearity 
   statistics, and casewise regression diagnostics.}
\usage{COUNT_REGRESSION(data, DV, forced = NULL, hierarchical = NULL,
                 family = 'poisson',
                 offset = NULL,
                 plot_type = 'residuals',
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

  \item{family}{
  \code{}(optional) The name of the error distribution to be 
  used in the model. The options are:
  \itemize{
   \item{"poisson" (the default),}
   \item{"quasipoisson", or}
   \item{"negbin", for negative binomial.}
  }
  Example: family = 'quasipoisson'}

  \item{offset}{
  \code{}(optional) The name of the offset variable, if there is one. This variable
  should be in the desired metric (e.g., log). No transformation of an
  offset variable is performed internally.
  \cr Example: offset = 'Varname'}

  \item{plot_type}{
  \code{}(optional) The kind of plots, if any. The options are:
  \itemize{
   \item{'residuals' (the default),}
   \item{'diagnostics', for regression diagnostics, and}
   \item{'none', for no plots.}
  }
  Example: plot_type = 'diagnostics'}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE. If TRUE, plots of residuals are also produced.}
}
\details{
	This function uses the glm function from the stats package, and the 
	negative.binomial function from the MASS package, and
	supplements the output with additional statistics and in formats that resembles
	SPSS and SAS output. The predictor variables can be numeric or factors.
	
	Predicted values for this model, for selected levels of the predictor variables,
	can be produced and plotted using the PLOT_MODEL funtion in this package.

  Good sources for interpreting count data regression residuals and diagnostics plots:
	\itemize{
  \item{\href{https://rpubs.com/benhorvath/glm_diagnostics}{rpubs.com/benhorvath}}
  \item{\href{https://library.virginia.edu/data/articles/understanding-deviance-residuals}{library.virginia.edu}}
  \item{\href{https://online.stat.psu.edu/stat462/node/207/}{online.stat.psu.edu}}
  }
}
\value{An object of class "COUNT_REGRESSION". The object is a list containing the
	following possible components:
\item{modelMAIN}{All of the glm function output for the regression model.}
\item{modelMAINsum}{All of the summary.glm function output for the regression model.}
\item{modeldata}{All of the predictor and outcome raw data that were used in the model,
	along with regression diagnostic statistics for each case.}
\item{collin_diags}{Collinearity diagnostic coefficients for models without interaction terms.}
\item{cormat}{The correlation matrix for the model coefficients.}
}
\references{
       {Atkins, D. C., & Gallop, R. J. (2007). Rethinking how family researchers 
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
\cr\cr {Hardin, J. W., & Hilbe, J. M. (2007). \emph{Generalized linear models 
        and extensions.} Stata Press.}
\cr\cr {Orme, J. G., & Combs-Orme, T. (2009). \emph{Multiple regression with discrete 
        dependent variables.} Oxford University Press.}
\cr\cr {Rindskopf, D. (2023). Generalized linear models. In H. Cooper, M. N. 
        Coutanche, L. M. McMullen, A. T. Panter, D. Rindskopf, & K. J. Sher (Eds.), 
	     \emph{APA handbook of research methods in psychology: Data analysis and 
	     research publication, } (2nd ed., pp. 201-218). American Psychological Association.}
}
\author{Brian P. O'Connor}
\examples{
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='OVRJOYED', 
                 forced=c('AGE','EDUC','REALRINC','SEX_factor'))
\donttest{
# negative binomial regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='HURTATWK', 
                 forced=c('AGE','EDUC','REALRINC','SEX_factor'),
                 family = 'negbin',
                 plot_type = 'diagnostics')

# with an offset variable
COUNT_REGRESSION(data=data_Orme_2009_5, DV='NumberAdopted', forced=c('Married'), 
                 offset='lnYearsFostered')
}
}
