\name{PARTIAL_COR}
\alias{PARTIAL_COR}
\title{Partial and semipartial correlations}
\description{Produces partial correlations between two or more variables (in set Y)
    while statistically controlling for one or more covariates (set C).
    It also produces partial correlations, semipartial correlations, and standardized 
    regression coefficients for predicting variables (in set Y) from one or 
    more set X variables.
}
\usage{PARTIAL_COR(data, Y, X=NULL, C=NULL, Ncases=NULL, verbose=TRUE)
}
\arguments{
  \item{data}{
  \code{}Either a dataframe of raw data (where the rows are cases and the columns are 
  the variables), or a square correlation matrix with row and column names.}
  	    
  \item{Y}{
  \code{}The names of one or more continuous variables in data.
  \cr Example: Y = c('var1', 'var2', 'var3')}

  \item{C}{
  \code{}The names of one or more continuous variables in data to be partialled
  out of the Y variable correlations. 
  \cr Example: C = c('var4', 'var5')}
  	      	                     
  \item{X}{
  \code{}The names of one or more continuous predictor variables in data. 
  \cr Example: X = c('var6', 'var7', 'var8')}
  
 \item{Ncases}{
  \code{}The number of cases. Required only when the input (data) is a correlation matrix.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? 
  \cr The options are: TRUE (default) or FALSE.}
}
\details{Y must be provided along with either one, or both, of C and X.


\strong{If Y and C are provided, but not X,} then the function computes:
  	\itemize{
  	\item{the correlations between the Y variables after partialling the 
  	C variables out of the Y variables.}
  }


  \strong{If Y and X are provided, but not C,} then the function computes:
  	\itemize{
  	\item{the standardized betas for the X variables predicting the Y variables;}
  
  	\item{the partial correlations for the X variables predicting the Y variables.
  In other words, for any given X variable, the other X variables are partialled
  out of both the given X variable and the Y variables. And,}
  
  	\item{the semi-partial correlations for the X variables predicting the Y variables.
  In other words, for any given X variable, the other X variables are partialled
  out of the given X variable and the Y variables remain as they are, untouched.}
  }


  \strong{If Y, X, and C are provided,} then the function computes:
  	\itemize{
  	\item{the correlations between the Y variables after partialling the C variables
  out of the Y variables;}

  	\item{the betas for the X variables predicting the C-partialled Y variables;}

  	\item{the partial correlations for the X variables predicting the C-partialled 
  Y variables. In other words, for any given X variable, the other X variables 
  are partialled out of both the given X variable and the C-partialled Y variables.
  And,}
  
  	\item{the semi-partial correlations for the X variables predicting the 
  C-partialled Y variables. In other words, for any given X variable, the 
  other X variables are partialled out of the given X variable but not out of the 
  C-partialled Y variables.}
  }
}
\value{A list containing the correlations, standardized regression coefficients (betas), 
	partial correlations, semi-partial correlations, t-test values, and p values.
}
\references{
      {Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	    multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	    Lawrence Erlbaum Associates.}
}
\author{Brian P. O'Connor }
\examples{
# C, but no X variables
PARTIAL_COR(data = data_DeLeo_2013, 
            Y = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
            C = c('Age','Parents_Income'), 
            X = NULL)

# X, but no C variables
PARTIAL_COR(data = data_DeLeo_2013, 
            Y = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
            C = NULL, 
            X = c('Impulsivity','Social_Interaction_Anxiety',
                  'Social_Support','Intolerance_of_Deviance','Family_Morals',
                  'Grade_Point_Average','Depression','Family_Conflict'))

# both X & C variables
PARTIAL_COR(data = data_DeLeo_2013, 
            Y = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use'), 
            C = c('Age','Parents_Income'), 
            X = c('Impulsivity','Social_Interaction_Anxiety',
                  'Social_Support','Intolerance_of_Deviance','Family_Morals',
                  'Grade_Point_Average','Depression','Family_Conflict'))
}
