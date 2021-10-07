\name{data_Lorah_Wong_2018}
\alias{data_Lorah_Wong_2018}
\docType{data}
\title{data_Lorah_Wong_2018}
\description{Moderated regression data from Lorah and Wong (2018). 
}
\usage{data(data_Lorah_Wong_2018)}

\source{
Lorah, J. A. & Wong, Y. J. (2018). Contemporary applications of moderation 
analysis in counseling psychology. \emph{Journal of Counseling Psychology, 65(5),} 629-640.
}

\examples{
head(data_Lorah_Wong_2018)

\donttest{
SIMPLE.REGRESSION(data=data_Lorah_Wong_2018, DV='sis',
       IV='burden', IV_type = 'numeric', IV_range='tumble',
       MOD='belong', MOD_type = 'numeric', MOD_levels='quantiles', MOD_range=NULL,
       quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
       CENTER = TRUE, COVARS='dep',
       PLOT_type = 'regions', PLOT_title=NULL,
       Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
       JN_type = 'Huitema', verbose=TRUE ) 
}
}
