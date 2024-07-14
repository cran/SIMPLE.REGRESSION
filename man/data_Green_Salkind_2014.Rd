\name{data_Green_Salkind_2014}
\alias{data_Green_Salkind_2014}
\docType{data}
\title{data_Green_Salkind_2014}
\description{Mutiple regression data from Green and Salkind (2018). 
}
\usage{data(data_Green_Salkind_2014)}
\source{
Green, S. B., & Salkind, N. J. (2014). Lesson 34: Multiple linear regression
(pp. 257-269). In, \emph{Using SPSS for Windows and Macintosh: Analyzing and 
understanding data.} New York, NY: Pearson.
}
\examples{
head(data_Green_Salkind_2014)
\donttest{
# forced (simultaneous) entry; replicating the output on p. 263	
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               forced=c('quads','gluts','abdoms','arms','grip')) 

# hierarchical entry; replicating the output on p. 265-266	
OLS_REGRESSION(data=data_Green_Salkind_2014, DV='injury', 
               hierarchical = list( step1=c('quads','gluts','abdoms'), 
                                    step2=c('arms','grip')) )
}
}
