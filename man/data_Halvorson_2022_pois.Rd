\name{data_Halvorson_2022_pois}
\alias{data_Halvorson_2022_pois}
\docType{data}
\title{data_Halvorson_2022_pois}
\description{Poisson regression data from Halvorson et al. (2022, p. 293). 
}
\usage{data(data_Halvorson_2022_pois)}
\source{
Halvorson, M. A., McCabe, C. J., Kim, D. S., Cao, X., & King, K. M. (2022). 
Making sense of some odd ratios: A tutorial and improvements to present practices 
in reporting and visualizing quantities of interest for binary and count outcome 
models. \emph{Psychology of Addictive Behaviors, 36(3),} 284-295.
}
\examples{
head(data_Halvorson_2022_pois)
\donttest{
# replicating Table 3, p 293
pois_Halvorson <-
  COUNT_REGRESSION(data=data_Halvorson_2022_pois, DV='Neg_OH_conseqs', 
                   forced=c('Gender','Positive_urgency','Planning',
                            'Sensation_seeking'), 
                   plot_type = 'diagnostics')

# replicating Figure 4, p 294
PLOT_MODEL(modobject = pois_Halvorson, 
           IV_focal_1 = 'Positive_urgency',   
           IV_focal_2 = 'Gender',
           bootstrap=FALSE, N_sims=1000, CI_level=95, 
           ylim = c(0, 20), 
           xlab = 'Positive Urgency',
           ylab = 'Expected Count of Alcohol Consequences', 
           title = 'Expected Count of Alcohol Consequences 
                    by Positive Urgency and Gender') 
}
}
