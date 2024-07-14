\name{data_Pituch_Stevens_2016}
\alias{data_Pituch_Stevens_2016}
\docType{data}
\title{data_Pituch_Stevens_2016}
\description{Logistic regression data from Pituch and Stevens (2016), Chapter 11. 
}
\usage{data(data_Pituch_Stevens_2016)}
\source{
Pituch, K. A., & Stevens, J. P. (2016). 
	      \emph{Applied multivariate statistics for the social sciences : Analyses with 
      	SAS and IBMs SPSS,} (6th ed.). Routledge.
}
\examples{
LOGISTIC_REGRESSION(data = data_Pituch_Stevens_2016, DV='Health', 
                    forced= c('Treatment','Motivation'))
}

