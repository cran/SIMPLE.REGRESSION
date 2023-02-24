\name{data_Chapman_Little_2016}
\alias{data_Chapman_Little_2016}
\docType{data}
\title{data_Chapman_Little_2016}

\description{Moderated regression data from Chapman and Little (2016).}

\usage{data(data_Chapman_Little_2016)}

\source{
Chapman, D. A., & Little, B. (2016). Climate change and disasters: How framing affects 
justifications for giving or withholding aid to disaster victims. 
\emph{Social Psychological and Personality Science, 7,} 13-20.
}

\examples{
head(data_Chapman_Little_2016)

\donttest{ 
# the data used by Hayes (2018, Introduction to Mediation, Moderation, and 
# Conditional Process Analysis: A Regression-Based Approach), replicating p. 239
SIMPLE.REGRESSION(data=data_Chapman_Little_2016, DV='justify',
       IV='frame', IV_type = 'numeric', IV_range='tumble',
       MOD='skeptic', MOD_type = 'numeric', MOD_levels='AikenWest', MOD_range=NULL,
       quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
       CENTER = FALSE, COVARS=NULL,
       PLOT_type = 'regions', PLOT_title=NULL, 
       Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
       JN_type = 'Huitema', verbose=TRUE ) 
}
}
