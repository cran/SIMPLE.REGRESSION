\name{data_Bauer_Curran_2005}
\alias{data_Bauer_Curran_2005}
\docType{data}
\title{data_Bauer_Curran_2005}
\description{Multilevel moderated regression data from Bauer and Curran (2005). 
}
\usage{data(data_Bauer_Curran_2005)}
\source{
Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel 
regression: Inferential and graphical techniques. \emph{Multivariate Behavioral 
Research, 40(3),} 373-400.
}
\examples{
head(data_Bauer_Curran_2005)
\donttest{
HSBmod <-nlme::lme(MathAch ~ Sector + CSES + CSES:Sector,
                   data = data_Bauer_Curran_2005, 
                   random = ~1 + CSES|School, method = "ML") 
summary(HSBmod)

REGIONS_OF_SIGNIFICANCE(model=HSBmod,  
                        plot_title='Johnson-Neyman Regions of Significance', 
                        Xaxis_label='Child SES',
                        Yaxis_label='Slopes of School Sector on Math achievement')  
}
}
