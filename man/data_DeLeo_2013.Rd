\name{data_DeLeo_2013}
\alias{data_DeLeo_2013}
\docType{data}
\title{data_DeLeo_2013}
\description{A dataset with multiple continuous variables that simulate the data
from De Leo and Wulfert (2013). The dataset is used in the examples for the present
PARTIAL_COR and SET_CORRELATION functions.
}
\usage{data(data_DeLeo_2013)}
\source{
	{De Leo, J. A., & Wulfert, E. (2013). Problematic internet use and other risky 
	behaviors in college students: An application of problem-behavior theory. 
	\emph{Psychology of Addictive Behaviors, 27(1),} 133-141.}
}
\examples{
head(data_DeLeo_2013)
\donttest{
# bipartial
SET_CORRELATION(data=data_DeLeo_2013, 
                IVs = c('Grade_Point_Average','Family_Morals','Social_Support',
                        'Intolerance_of_Deviance','Impulsivity',
                        'Social_Interaction_Anxiety'), 
                DVs = c('Problematic_Internet_Use','Tobacco_Use','Alcohol_Use',
                        'Illicit_Drug_Use'), 
                IV_covars = c('Age','Parents_Income'), 
                DV_covars = c('Gambling_Behavior','Unprotected_Sex'),
                display_cormats=TRUE) }
}
