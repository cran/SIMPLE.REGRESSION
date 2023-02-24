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
	Lawrence Erlbaum Associates Publishers.
}

\examples{
head(data_Cohen_Aiken_West_2003_7)

\donttest{
# replicates p 276 of Chapter 7 of Cohen, Cohen, West, & Aiken (2003)
SIMPLE.REGRESSION(data=data_Cohen_Aiken_West_2003_7, DV='yendu',
                  IV='xage', IV_type = 'numeric', IV_range='tumble',
                  MOD='zexer', MOD_type = 'numeric', MOD_levels='AikenWest', 
                  MOD_range=NULL,
                  quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                  CENTER = TRUE, COVARS=NULL, 
                  PLOT_type = 'regions', PLOT_title=NULL, 
                  Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
                  JN_type = 'Huitema', verbose=TRUE ) 
}
}
