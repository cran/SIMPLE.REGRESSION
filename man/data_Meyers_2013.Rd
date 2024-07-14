\name{data_Meyers_2013}
\alias{data_Meyers_2013}
\docType{data}
\title{data_Meyers_2013}
\description{Logistic regression data from Myers et al. (2013). 
}
\usage{data(data_Meyers_2013)}
\source{
Meyers, L. S., Gamst, G. C., & Guarino, A. J. (2013). Chapter 30: Binary 
logistic regression. \emph{Performing data analysis using IBM SPSS.}
Hoboken, NJ: Wiley.
}
\examples{
head(data_Meyers_2013)
\donttest{
LOGISTIC_REGRESSION(data= data_Meyers_2013, DV='graduated', forced= c('sex','family_encouragement'))
}
}
