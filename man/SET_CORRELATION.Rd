\name{SET_CORRELATION}
\alias{SET_CORRELATION}
\title{Cohen's Set Correlation Analysis}
\description{Performs Cohen's set correlation analysis of associations
between two sets of variables while statistically controlling for one or more
other variables. Estimates of overall, multivariate association between the
two sets of variables are provided, along with partial correlations and output
from OLS regression analyses for each dependent variable.}
\usage{SET_CORRELATION(data, IVs, DVs, IV_covars=NULL, DV_covars=NULL,
                Ncases=NULL, verbose=TRUE, display_cormats=FALSE)
}
\arguments{
  \item{data}{
  \code{}Either a dataframe of raw data (where the rows are cases and the columns are 
  the variables), or a square correlation matrix with row and column names.}

  \item{IVs}{
  \code{}The name(s) of the independent/predictor variable(s) in data. 
  \cr Example: IVs = c('var1', 'var2', 'var3')}
  
  \item{DVs}{
  \code{}The name(s) of the dependent variable(s) in data. 
  \cr Example: DVs = c('var4', 'var5', 'var6')}

  \item{IV_covars}{
  \code{}The name(s) of the variable(s), if any, to be partialled out of the IVs. 
  \cr Example: IV_covars = c('var7', 'var8')}

  \item{DV_covars}{
  \code{}The name(s) of the variable(s), if any, to be partialled out of the DVs. 
  \cr Example: DV_covars = c('var9', 'var10')}

 \item{Ncases}{
  \code{}The number of cases. Required only when the input (data) is a correlation matrix.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? The options are: 
  TRUE (default) or FALSE.}

  \item{display_cormats}{
  \code{}Should the variable correlation matrices be displayed in console? 
  The options are: TRUE or FALSE(default).}
}
\details{
	Set correlation analysis and canonical correlation analysis are both fully 
	multivariate methods for examining associations between two sets of variables.
	However, in CCA the focus is on linear combinations of predictor and
	criterion variables, which are often difficult to interpret.
	In contrast, in set correlation analysis the focus is typically on the
	associations between two sets of variables while statistically
	controlling for other variables (rather than on linear combinations).
	The outcome variables of interest in set correlation analysis are the
	(possibly partialled) dependent variables themselves and not 
	composites of variables.
	
	A key feature of set correlation analysis is the option of examining
	the overlap between two sets of variables while statistically controlling 
	for one or more other variables. The covariates that are removed from one
	set of variables (e.g., the DVs) may or may not be the same covariates
	that are removed from the other set of variables (e.g., the IVs).
	
	In the present function, when there is a wish to statistically remove
	the same covariates from both sets (i.e., from both the IVs and DVs), then 
	simply enter the same covariate names on both the IV_covars and DV_covars
	arguments.
	
	The options together result in five different types of data scenarios that
	can be examined:
	
	\strong{Whole}, in which the associations between two sets (IVs and DVs) are assessed
	without any partialling out whatsoever;
	
	\strong{Partial}, in which the associations between two sets (IVs and DVs) are assessed
	while partialling the same covariates (one or more) out of both the IVs and DVs;
	
	\strong{X Semipartial}, in which the associations between two sets (IVs and DVs) are assessed
	while partialling one or more covariates out of the IV set while leaving the variables
	in the DV set untouched (unpartialled);
	
	\strong{Y Semipartial}, in which the associations between two sets (IVs and DVs) are assessed
	while partialling one or more covariates out of the DV set while leaving the variables
	in the IV set untouched (unpartialled); and
	
	\strong{Bipartial}, in which the associations between two sets (IVs and DVs) are assessed
	while partialling one or more covariates out of the DV set and while partialling 
	one or more other (different) covariates out of the IV set.
	
	The set correlation analyses in this function are conducted using only the
	correlations between the variables. When raw data are entered into the function,
	the variable correlation matrix is computed and becomes the sole basis of
	all further set correlation analyses.
}
\value{An object of class "SET_CORRELATION". The object is a list containing the
	following components:
\item{bigR}{The Pearson correlation matrix for the variables in the analyses.}
\item{Ryy}{The correlations between the DVs.}
\item{Rxx}{The correlations between the IVs.}
\item{Rx_y}{The correlation between the DVs and IVs.}
\item{betas}{The standardized betas.}
\item{se_betas}{The standard errors of the standardized betas.}
\item{t}{The t test values for the standardized betas.}
\item{pt}{The p values for the t tests for the standardized betas.}
}
\references{
       {Cohen, J. (1982). Set correlation as a general multivariate data-analytic 
       method. \emph{Multivariate Behavioral Research, 17(3)}, 301-341.}
\cr\cr {Cohen, J. (1988). Set correlation and multivariate Methods.
       In J. Cohen, \emph{Statistical power analysis for the behavioral 
       sciences} (2nd ed., pp. 467-530). Mahwah, NJ: Erlbaum.}
\cr\cr {Cohen, J. (1993). Set correlation. In G. Keren & C. Lewis (Eds.), \emph{A 
        handbook for data analysis in the behavioral sciences: Statistical 
        issues} (pp. 165-198). Mahwah, NJ: Erlbaum.}
\cr\cr {Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). Multiple 
        dependent variables: Set correlation. In, \emph{Applied 
	      multiple regression/correlation analysis for the behavioral sciences} 
	      (3rd ed., pp. 608-628). Lawrence Erlbaum Associates.}
}
\author{Brian P. O'Connor}
\examples{
# data from Cohen et al. (2003)
Cohen_2003_p621 <- '
 1.0
 .53  1.0  
 .62  .61  1.0
 .19  .23  .03  1.0
-.09  .10  .10 -.02  1.0
 .08  .18  .12  .02  .05  1.0
 .02  .02  .03  .00  .06  .22  1.0
-.12 -.10 -.06 -.02  .18 -.07 -.01  1.0
 .08  .15  .12 -.02  .02  .36 -.05 -.03  1.0'

Cohen_2003_p621_noms <- c('ADHD', 'CD', 'ODD', 'Sex', 'Age', 'MONLY', 
                          'MWORK', 'MAGE', 'Poverty')

Cohen_2003_p621 <- data.matrix( read.table(text=Cohen_2003_p621, fill=TRUE, 
                                           col.names=Cohen_2003_p621_noms,
                                           row.names=Cohen_2003_p621_noms ))
Cohen_2003_p621[upper.tri(Cohen_2003_p621)] <- 
  t(Cohen_2003_p621)[upper.tri(Cohen_2003_p621)]

# whole
SET_CORRELATION(data=Cohen_2003_p621, 
                IVs = c('Sex', 'Age', 'MONLY', 'MWORK', 'MAGE', 'Poverty'), 
                DVs = c('ADHD', 'CD', 'ODD'), 
                IV_covars = NULL, 
                DV_covars = NULL,
                Ncases = 701) 
\donttest{
# bipartial
SET_CORRELATION(data=data_DeLeo_2013, 
                IVs = c('Grade_Point_Average','Family_Morals','Social_Support',
                        'Intolerance_of_Deviance','Impulsivity','Social_Interaction_Anxiety'), 
                DVs = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
                IV_covars = c('Age','Parents_Income'), 
                DV_covars = c('Gambling_Behavior','Unprotected_Sex'),
                display_cormats=TRUE) 

# X semipartial
SET_CORRELATION(data=data_DeLeo_2013, 
                IVs = c('Grade_Point_Average','Family_Morals','Social_Support',
                        'Intolerance_of_Deviance','Impulsivity','Social_Interaction_Anxiety'), 
                DVs = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
                IV_covars = c('Age','Parents_Income'), 
                DV_covars = NULL) 

# partial
SET_CORRELATION(data=data_DeLeo_2013, 
                IVs = c('Grade_Point_Average','Family_Morals','Social_Support',
                        'Intolerance_of_Deviance','Impulsivity','Social_Interaction_Anxiety'), 
                DVs = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
                IV_covars = c('Age','Parents_Income'), 
                DV_covars = c('Age','Parents_Income')) 
}
}
