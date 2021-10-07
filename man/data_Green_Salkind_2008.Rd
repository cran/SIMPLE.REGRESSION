\name{data_Green_Salkind_2008}
\alias{data_Green_Salkind_2008}
\docType{data}
\title{data_Green_Salkind_2008}
\description{Mutiple regression data from Green and Salkind (2018). 
}
\usage{data(data_Green_Salkind_2008)}

\source{
Green, S. B., & Salkind, N. J. (2008). Lesson 34: Multiple linear regression
(pp. 285-299). In, \emph{Using SPSS for Windows and Macintosh: Analyzing and 
understanding data.} New York, NY: Pearson.
}

\examples{
head(data_Green_Salkind_2008)

\donttest{
SIMPLE.REGRESSION(data=data_Green_Salkind_2008, DV='injury', 
                  forced=c('quads','gluts','abdoms','arms','grip')) 
}
}
