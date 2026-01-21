\name{data_Orme_2009_5}
\alias{data_Orme_2009_5}
\docType{data}
\title{data_Orme_2009_5}
\description{Data for count regression from Orme and Combs-Orme (2009), Chapter 5. 
}
\usage{data(data_Orme_2009_5)}
\source{
Orme, J. G., & Combs-Orme, T. (2009). \emph{Multiple Regression With Discrete 
        Dependent Variables.} Oxford University Press.
}
\examples{
head(data_Orme_2009_5)
\donttest{
# for Orme, 2009, p. 175: negative binomial regression with an offset variable
COUNT_REGRESSION(data=data_Orme_2009_5, DV='NumberAdopted', 
                 forced=c('Married','zParentRole'), 
                 model_type = 'negbin', 
                 offset='lnYearsFostered')
}
}

