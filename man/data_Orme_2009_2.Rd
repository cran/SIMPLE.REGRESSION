\name{data_Orme_2009_2}
\alias{data_Orme_2009_2}
\docType{data}
\title{data_Orme_2009_2}
\description{Logistic regression data from Orme and Combs-Orme (2009), Chapter 2. 
}
\usage{data(data_Orme_2009_2)}
\source{
Orme, J. G., & Combs-Orme, T. (2009). \emph{Multiple Regression With Discrete 
        Dependent Variables.} Oxford University Press.
}
\examples{
head(data_Orme_2009_2)
\donttest{
# for Orme, 2009, p. 60
LOGISTIC_REGRESSION(data = data_Orme_2009_2, DV='ContinueFostering', 
                    forced= c('zResources', 'Married'))
}
}

