\name{data_Cohen_Aiken_West_2003_9}
\alias{data_Cohen_Aiken_West_2003_9}
\docType{data}
\title{data_Cohen_Aiken_West_2003_9}
\description{Moderated regression data for a continuous predictor and a
	categorical moderator from Cohen, Cohen, West, & Aiken (2003, Chapter 9). 
}
\usage{data(data_Cohen_Aiken_West_2003_9)}
\source{
Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	Lawrence Erlbaum Associates.
}
\examples{
head(data_Cohen_Aiken_West_2003_9)
\donttest{
# replicates p 376 of Chapter 9 of Cohen, Cohen, West, & Aiken (2003)
MODERATED_REGRESSION(data=data_Cohen_Aiken_West_2003_9, DV='SALARY',
                     IV='PUB', IV_range='tumble',
                     MOD='DEPART_f', MOD_type = 'factor', MOD_levels='AikenWest', 
                     MOD_reflevel = 'psychology',
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = TRUE,  
                     plot_type = 'interaction') 
}
}
