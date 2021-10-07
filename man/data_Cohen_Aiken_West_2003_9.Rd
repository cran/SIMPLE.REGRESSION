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
	Lawrence Erlbaum Associates Publishers.
}

\examples{
head(data_Cohen_Aiken_West_2003_9)

\donttest{
SIMPLE.REGRESSION(data=data_Cohen_Aiken_West_2003_9, DV='SALARY', forced=c('PUB','DEPART_f'))

# replicates p 376 of Chapter 9 of Cohen, Cohen, West, & Aiken (2003)
SIMPLE.REGRESSION(data=data_Cohen_Aiken_West_2003_9, DV='SALARY', forced=NULL,
                  IV='PUB', IV_type = 'numeric', IV_range='tumble',
                  MOD='DEPART_f', MOD_type = 'factor', MOD_levels='AikenWest', 
                  MOD_range=NULL,
                  quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                  CENTER = TRUE, COVARS=NULL, 
                  PLOT_type = 'regions', PLOT_title=NULL, 
                  Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
                  JN_type = 'Huitema', verbose=TRUE ) 
}
}
