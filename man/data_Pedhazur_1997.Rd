\name{data_Pedhazur_1997}
\alias{data_Pedhazur_1997}
\docType{data}
\title{data_Pedhazur_1997}
\description{Moderated regression data for a continuous predictor and a
	dichotomous moderator from Pedhazur (1997, p. 588). 
}
\usage{data(data_Pedhazur_1997)}
\source{
Pedhazur, E. J. (1997). \emph{Multiple regression in behavioral research: Explanation 
   and prediction.} (3rd ed.). Fort Worth, Texas: Wadsworth Thomson Learning.
}
\examples{
head(data_Pedhazur_1997)
\donttest{
# replicating results on p. 594 for the Johnson-Neyman technique for a categorical moderator	
MODERATED_REGRESSION(data=data_Pedhazur_1997, DV='Y', 
                     IV='X', IV_range='tumble',
                     MOD='Directive', MOD_type = 'factor', MOD_levels='quantiles', 
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = FALSE, 
                     plot_type = 'interaction', 
                     JN_type = 'Pedhazur') 
}
}
