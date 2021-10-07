\name{PARTIAL_COEFS}
\alias{PARTIAL_COEFS}
\title{Standardized coefficients and partial correlations for multiple regression}
\description{Produces standardized regression coefficients, partial correlations,   
	and semi-partial correlations for a correlation matrix in which one variable 
    is a dependent or outcome variable and the other variables are independent
    or predictor variables.}
\usage{PARTIAL_COEFS(cormat, modelRsq=NULL, verbose=TRUE)}
\arguments{
  \item{cormat}{
  \code{}A correlation matrix. The DV (the dependent or outcome variable) must be in 
  the first row/column of cormat.}
  	    
  \item{modelRsq}{
  \code{}optional. The model Rsquared, which makes the computations slightly faster
  when it is available.}
  	      	                     
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{An object of class "data.frame" containing with the standardized regression coefficients (betas), 
	the Pearson correlations, the partial correlations, and the semi-partial correlations for each
	variable with the DV.
}
\references{
       {Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	    multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	    Lawrence Erlbaum Associates Publishers.}
}
\author{Brian P. O'Connor }
\examples{
PARTIAL_COEFS(cormat = cor(data_Green_Salkind_2008))
}
