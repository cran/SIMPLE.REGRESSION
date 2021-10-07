\name{REGIONS_OF_SIGNIFICANCE}
\alias{REGIONS_OF_SIGNIFICANCE}
\title{Plots of Johnson-Neyman regions of significance for interactions}

\description{Plots of Johnson-Neyman regions of significance for interactions 
	 in moderated multiple regression, for both SIMPLE.REGRESSION models (objects) and 
	 for lme models from the nlme package.}
	
\usage{REGIONS_OF_SIGNIFICANCE(model,
                        IV_range=NULL, MOD_range=NULL,
                        PLOT_title=NULL, Xaxis_label=NULL,
                        Yaxis_label=NULL, LEGEND_label=NULL,
                        namesIVMOD_raw=NULL, namesIVMOD_model=NULL)
}
\arguments{
  \item{model}{
  \code{}The name of a SIMPLE.REGRESSION model, or of an lme model from the nlme package.}

  \item{IV_range}{
  \code{}(optional) The range of the IV to be used in the plot.}
  
  \item{MOD_range}{
  \code{}(optional)  The range of the MOD values to be used in the }
  
  \item{PLOT_title}{
  \code{}(optional)  The plot title.}
  
  \item{Xaxis_label}{
  \code{}(optional)  A label for the X axis to be used in the requested plot.}
  
  \item{Yaxis_label}{
  \code{}(optional)  A label for the Y axis to be used in the requested plot.}
  
  \item{LEGEND_label}{
  \code{}(optional)  The legend label for a moderated regression.}
    
  \item{namesIVMOD_raw}{
  \code{}optional, and for lme/nlme models only. If model is an lme object & IV 
  is a two-level factor, then namesIVMOD_model must be specified (because lme alters 
  the variable names).}
  
  \item{namesIVMOD_model}{
  \code{}optional, and for lme/nlme models only. The namesIVMOD_model argument can 
  be used to id the key terms () from an lme model that involved more than IV, 
  MOD, & Xn terms. The argument is used only to create the key B and S objects 
  for the J-N analyses. Other terms in the model are ignored.}
  
  % \item{verbose}{
  % \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) 
           % or FALSE.}
}
\value{An object of class "SIMPLE.REGRESSION". The object is a list containing the
	following possible components:
\item{JN.data}{The Johnson-Neyman results for a moderated regression.}
\item{ros}{The Johnson-Neyman regions of significance for a moderated regression.}
}

\references{
    {Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel 
	regression: Inferential and graphical techniques. \emph{Multivariate Behavioral 
	Research, 40(3),} 373-400.}
\cr\cr {Huitema, B. (2011). \emph{The analysis of covariance and alternatives: Statistical  
    methods for experiments, quasi-experiments, and single-case studies.} Hoboken, NJ: Wiley.}
\cr\cr {Johnson, P. O., & Neyman, J. (1936). Tests of certain linear hypotheses and their 
    application to some educational problems. \emph{Statistical Research Memoirs, 1,} 57-93.}
\cr\cr {Johnson, P. O., & Fey, L. C. (1950). The Johnson-Neyman technique, its theory, and 
	application. \emph{Psychometrika, 15,} 349-367.}
\cr\cr {Pedhazur, E. J. (1997). \emph{Multiple regression in behavioral research: Explanation 
    and prediction.} (3rd ed.). Fort Worth, Texas: Wadsworth Thomson Learning}
\cr\cr {Rast, P., Rush, J., Piccinin, A. M., & Hofer, S. M. (2014). The identification of regions of 
	significance in the effect of multimorbidity on depressive symptoms using longitudinal data: An 
	application of the Johnson-Neyman technique. \emph{Gerontology, 60,} 274-281.}
}
\author{Brian P. O'Connor }
\examples{
head(data_Cohen_Aiken_West_2003_7)

CAW_7 <- 
SIMPLE.REGRESSION(data=data_Cohen_Aiken_West_2003_7, DV='yendu',
                  IV='xage', IV_type = 'numeric', IV_range='tumble',
                  MOD='zexer', MOD_type = 'numeric', MOD_levels='quantiles', 
                  MOD_range=NULL,
                  quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                  CENTER = TRUE, COVARS=NULL, 
                  PLOT_type = 'interaction', PLOT_title=NULL, 
                  Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
                  JN_type = 'Huitema', verbose=TRUE ) 

REGIONS_OF_SIGNIFICANCE(model=CAW_7,  
        IV_range=NULL, MOD_range='minmax', 
        PLOT_title=NULL, Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
        namesIVMOD_raw=NULL, namesIVMOD_model=NULL) 
	

\donttest{
head(data_Bauer_Curran_2005)

HSBmod <- nlme::lme(MathAch ~ Sector + CSES + CSES:Sector , 
			        data = data_Bauer_Curran_2005, 
			        random = ~1 + CSES|School, method = "ML") 
summary(HSBmod)

REGIONS_OF_SIGNIFICANCE(model=HSBmod,  
        IV_range=NULL, MOD_range=NULL, 
        PLOT_title=NULL, Xaxis_label=NULL, Yaxis_label=NULL, LEGEND_label=NULL,
        namesIVMOD_raw=NULL, namesIVMOD_model=NULL) 
}
}
