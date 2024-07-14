\name{data_Huitema_2011}
\alias{data_Huitema_2011}
\docType{data}
\title{data_Huitema_2011}
\description{Moderated regression data for a continuous predictor and a
	dichotomous moderator from Huitema (2011, p. 253). 
}
\usage{data(data_Huitema_2011)}
\source{
Huitema, B. (2011). \emph{The analysis of covariance and alternatives: Statistical methods 
for experiments, quasi-experiments, and single-case studies.} Hoboken, NJ: Wiley.
}
\examples{
head(data_Huitema_2011)
\donttest{
# replicating results on p. 255 for the Johnson-Neyman technique for a categorical moderator
MODERATED_REGRESSION(data=data_Huitema_2011, DV='Y', 
                     IV='X', IV_range='tumble',
                     MOD='D', MOD_type = 'factor',  
                     center = FALSE,  
                     plot_type = 'interaction',
                     JN_type = 'Huitema') 
}
}
