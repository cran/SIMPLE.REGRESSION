\name{data_Cohen_Aiken_West_2003_7}
\alias{data_Cohen_Aiken_West_2003_7}
\docType{data}
\title{data_Cohen_Aiken_West_2003_7}
\description{Moderated regression data for a continuous predictor and a
	continuous moderator from Cohen, Cohen, West, & Aiken (2003, Chapter 7). 
}
\usage{data(data_Cohen_Aiken_West_2003_7)}
\source{
Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied 
	multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). 
	Lawrence Erlbaum Associates.
}
\examples{
head(data_Cohen_Aiken_West_2003_7)
\donttest{
# replicates p 276 of Chapter 7 of Cohen, Cohen, West, & Aiken (2003)
MODERATED_REGRESSION(data=data_Cohen_Aiken_West_2003_7, DV='yendu',
                     IV='xage', IV_range='tumble',
                     MOD='zexer', MOD_levels='AikenWest', 
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     center = TRUE, 
                     plot_type = 'regions') 
}
}
