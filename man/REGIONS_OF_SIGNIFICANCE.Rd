\name{REGIONS_OF_SIGNIFICANCE}
\alias{REGIONS_OF_SIGNIFICANCE}
\title{Plots of Johnson-Neyman regions of significance for interactions}
\description{Plots of Johnson-Neyman regions of significance for interactions 
	 in moderated multiple regression, for both MODERATED_REGRESSION models (which are
	 produced by this package) and for lme models (from the nlme package).
}
\usage{REGIONS_OF_SIGNIFICANCE(model,
                        IV_range=NULL, MOD_range=NULL,
                        plot_title=NULL, Xaxis_label=NULL,
                        Yaxis_label=NULL, legend_label=NULL,
                        names_IV_MOD=NULL)
}
\arguments{
  \item{model}{
  \code{}The name of a MODERATED_REGRESSION model, or of an lme model from the nlme package.}

  \item{IV_range}{
  \code{}(optional) The range of the IV to be used in the plot.
  \cr Example: IV_range = c(1, 10)}
  
  \item{MOD_range}{
  \code{}(optional)  The range of the MOD values to be used in the plot.
  \cr Example: MOD_range = c(2, 4, 6)}
  
  \item{plot_title}{
  \code{}(optional)  The plot title.
  \cr Example: plot_title = 'Regions of Significance Plot'}
  
  \item{Xaxis_label}{
  \code{}(optional)  A label for the X axis to be used in the plot.
  \cr Example: Xaxis_label = 'IV name'}
  
  \item{Yaxis_label}{
  \code{}(optional)  A label for the Y axis to be used in the plot.
  \cr Example: Yaxis_label = 'DV name'}
  
  \item{legend_label}{
  \code{}(optional)  The legend label.
  \cr Example: legend_label = 'Simple Slopes'}
    
  \item{names_IV_MOD}{
  \code{}(optional) and for lme/nlme models only. Use this argument to ensure
  that the IV and MOD variables are correctly identified for the plot.
  There are three scenarios in particular that may require specification of this
  argument:
  \itemize{
    \item{when there are covariates in addition to IV & MOD as predictors,}
    \item{if the order of the variables in model is not IV then MOD, or},
    \item{if the IV is a two-level factor (because lme alters the variable names in this case).}
  }
  Example: names_IV_MOD = c('IV name', 'MOD name')}
}
\value{A list with the following possible components:
\item{JN.data}{The Johnson-Neyman results for a moderated regression.}
\item{ros}{The Johnson-Neyman regions of significance for a moderated regression.}
}
\references{
    {Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel 
	regression: Inferential and graphical techniques. \emph{Multivariate Behavioral 
	Research, 40(3),} 373-400.}
\cr\cr {Huitema, B. (2011). \emph{The analysis of covariance and alternatives: Statistical  
    methods for experiments, quasi-experiments, and single-case studies.} John Wiley & Sons.}
\cr\cr {Johnson, P. O., & Neyman, J. (1936). Tests of certain linear hypotheses and their 
    application to some educational problems. \emph{Statistical Research Memoirs, 1,} 57-93.}
\cr\cr {Johnson, P. O., & Fey, L. C. (1950). The Johnson-Neyman technique, its theory, and 
	application. \emph{Psychometrika, 15,} 349-367.}
\cr\cr {Pedhazur, E. J. (1997). \emph{Multiple regression in behavioral research: Explanation 
    and prediction.} (3rd ed.). Wadsworth Thomson Learning}
\cr\cr {Rast, P., Rush, J., Piccinin, A. M., & Hofer, S. M. (2014). The identification of regions of 
	significance in the effect of multimorbidity on depressive symptoms using longitudinal data: An 
	application of the Johnson-Neyman technique. \emph{Gerontology, 60,} 274-281.}
}
\author{Brian P. O'Connor }
\examples{
head(data_Cohen_Aiken_West_2003_7)

CAW_7 <- 
MODERATED_REGRESSION(data=data_Cohen_Aiken_West_2003_7, DV='yendu',
                     IV='xage',IV_range='tumble',
                     MOD='zexer', MOD_levels='quantiles', 
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     plot_type = 'interaction') 

REGIONS_OF_SIGNIFICANCE(model=CAW_7) 
\donttest{
head(data_Bauer_Curran_2005)

HSBmod <-nlme::lme(MathAch ~ Sector + CSES + CSES:Sector,
                   data = data_Bauer_Curran_2005, 
                   random = ~1 + CSES|School, method = "ML") 
summary(HSBmod)

REGIONS_OF_SIGNIFICANCE(model=HSBmod,  
                        plot_title='Johnson-Neyman Regions of Significance', 
                        Xaxis_label='Child SES',
                        Yaxis_label='Slopes of School Sector on Math achievement')  
                        

# moderated regression	-- with numeric values for IV_range & MOD_levels='AikenWest'       
mharsh_agg <- 
  MODERATED_REGRESSION(data=data_OConnor_Dvorak_2001, DV='Aggressive_Behavior',
                       IV='Maternal_Harshness', IV_range=c(1,7.7), 
                       MOD='Resiliency', MOD_levels='AikenWest', 
                       quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                       center = FALSE, 
                       plot_type = 'interaction', 
                       DV_range = c(1,6), 
                       Xaxis_label='Maternal Harshness', 
                       Yaxis_label='Adolescent Aggressive Behavior', 
                       legend_label='Resiliency') 

REGIONS_OF_SIGNIFICANCE(model=mharsh_agg,  
                        plot_title='Johnson-Neyman Regions of Significance', 
                        Xaxis_label='Resiliency', 
                        Yaxis_label='Slopes of Maternal Harshness on Aggressive Behavior') 
}
}
