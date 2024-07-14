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
model_Lorah <- 
MODERATED_REGRESSION(data=data_Lorah_Wong_2018, DV='suicidal',
                     IV='burden', IV_range='tumble',
                     MOD='belong_thwarted', MOD_levels='quantiles',
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     COVARS='depression', center = TRUE, 
                     plot_type = 'regions') 
       
REGIONS_OF_SIGNIFICANCE(model=model_Lorah,  
                        plot_title='Johnson-Neyman Regions of Significance', 
                        Xaxis_label='Thwarted Belongingness', 
                        Yaxis_label='Slopes of Burdensomeness on Suicical Ideation', 
                        legend_label=NULL)        
}
}
