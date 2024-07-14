\name{data_Kremelburg_2011}
\alias{data_Kremelburg_2011}
\docType{data}
\title{data_Kremelburg_2011}
\description{Logistic and Poisson regression data from Kremelburg (2011). 
}
\usage{data(data_Kremelburg_2011)}
\source{
Kremelburg, D. (2011). Chapter 6: Logistic, ordered, multinomial, negative
binomial, and Poisson regression. \emph{Practical statistics: A quick and easy guide 
to IBM SPSS Statistics, STATA, and other statistical software.} Sage.
}
\examples{
head(data_Kremelburg_2011)
\donttest{
LOGISTIC_REGRESSION(data = data_Kremelburg_2011, DV='OCCTRAIN',
                    hierarchical=list( step1=c('AGE'), step2=c('EDUC','REALRINC')) )
         
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='OVRJOYED', 
                 forced=c('AGE','EDUC','REALRINC','SEX_factor'))
}
}
