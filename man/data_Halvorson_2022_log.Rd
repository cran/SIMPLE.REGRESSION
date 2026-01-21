\name{data_Halvorson_2022_log}
\alias{data_Halvorson_2022_log}
\docType{data}
\title{data_Halvorson_2022_log}
\description{Logistic regression data from Halvorson et al. (2022, p. 291). 
}
\usage{data(data_Halvorson_2022_log)}
\source{
Halvorson, M. A., McCabe, C. J., Kim, D. S., Cao, X., & King, K. M. (2022). 
Making sense of some odd ratios: A tutorial and improvements to present practices 
in reporting and visualizing quantities of interest for binary and count outcome 
models. \emph{Psychology of Addictive Behaviors, 36(3),} 284-295.
}
\examples{
head(data_Halvorson_2022_log)
\donttest{
# replicating Figure 3, p 292
log_Halvorson <-
  LOGISTIC_REGRESSION(data=data_Halvorson_2022_log, DV='Y', forced=c('x1','x2'), 
                      plot_type = 'diagnostics')

# high & low values for x2
x2_high <- mean(data_Halvorson_2022_log$x1) + sd(data_Halvorson_2022_log$x1)
x2_low  <- mean(data_Halvorson_2022_log$x1) - sd(data_Halvorson_2022_log$x1)

PLOT_MODEL(modobject = log_Halvorson, 
           IV_focal_1 = 'x1',   
           IV_focal_2 = 'x2',  IV_focal_2_values = c(x2_low, x2_high),
           bootstrap=TRUE, N_sims=1000, CI_level=95, 
           ylim = c(0, 1), 
           xlab = 'x1',
           ylab = 'Expected Probability', 
           title = 'Probability of Y by x1 and x2 for Simulated Data Example') 
 
}
}
