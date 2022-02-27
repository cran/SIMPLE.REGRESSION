\name{data_OConnor_Dvorak_2001}
\alias{data_OConnor_Dvorak_2001}
\docType{data}
\title{data_OConnor_Dvorak_2001}
\description{Moderated regression data from O'Connor and Dvorak (2001)}

\details{A data frame with scores for 131 male adolescents on resiliency,
	maternal harshness, and aggressive behavior. The data are from
	O'Connor and Dvorak (2001, p. 17) and are provided as trial moderated 
	regression data for the SIMPLE_REGRESSION and REGIONS_OF_SIGNIFICANCE functions.
}

\references{
O'Connor, B. P., & Dvorak, T. (2001). Conditional associations between parental behavior 
and adolescent problems: A search for personality-environment interactions. 
\emph{Journal of Research in Personality, 35,} 1-26.
}

\examples{
head(data_OConnor_Dvorak_2001)

\donttest{
mharsh_agg <- SIMPLE.REGRESSION(data=data_OConnor_Dvorak_2001, DV='Aggressive_Behavior',
       IV='Maternal_Harshness', IV_type = 'numeric', IV_range=c(1,7.7), 
       MOD='Resiliency', MOD_type = 'numeric', MOD_levels='AikenWest', MOD_range=NULL,
       quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
       CENTER = FALSE, COVARS=NULL, 
       PLOT_type = 'interaction', 
       PLOT_title=NULL, 
       DV_range = c(1,6), 
       Xaxis_label='Maternal Harshness', 
       Yaxis_label='Adolescent Aggressive Behavior', 
       LEGEND_label='Resiliency',
       JN_type = 'Huitema', verbose=TRUE ) 

REGIONS_OF_SIGNIFICANCE(model=mharsh_agg,  
        IV_range=NULL, 
        MOD_range='minmax', 
        PLOT_title='Slopes of Maternal Harshness on Aggression by Resiliency', 
        Xaxis_label='Resiliency', 
        Yaxis_label='Slopes of Maternal Harshness on Aggressive Behavior ', 
        LEGEND_label=NULL,
        namesIVMOD_raw=NULL, 
        namesIVMOD_model=NULL) 
}
}
