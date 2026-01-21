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
# for Kremelburg, 2011, p. 244
LOGISTIC_REGRESSION(data = data_Kremelburg_2011, DV='OCCTRAIN',
                    forced = c('AGE','female','EDUC','REALRINC'))
                    
# for Kremelburg, 2011, p.262 & p. 282: poisson regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='OVRJOYED', 
                 forced=c('AGE','EDUC','REALRINC','female'))
                 
# for Kremelburg, 2011, p. 266-267 & p. 284: negative binomial regression
COUNT_REGRESSION(data=data_Kremelburg_2011, DV='HURTATWK', 
                 forced=c('AGE','EDUC','REALRINC','female'),
                 model_type = 'negbin')
}
}
