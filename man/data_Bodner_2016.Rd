\name{data_Bodner_2016}
\alias{data_Bodner_2016}
\docType{data}
\title{data_Bodner_2016}
\description{Moderated regression data used by Bodner (2016) to illustrate the
	tumble graphs method of plotting interactions. The data were also used by
	Bauer and Curran (2005). 
}
\usage{data(data_Bodner_2016)}
\source{
{Bodner, T. E. (2016). Tumble Graphs: Avoiding misleading end point 
	extrapolation when graphing interactions from a moderated multiple 
	regression analysis. 
	\emph{Journal of Educational and Behavioral Statistics, 41(6),} 593-604.}
\cr\cr {Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and 
	multilevel regression: Inferential and graphical techniques. 
	\emph{Multivariate Behavioral Research, 40(3),} 373-400.}
}
\examples{
head(data_Bodner_2016)
\donttest{
# replicates p 599 of Bodner (2016)
MODERATED_REGRESSION(data=data_Bodner_2016, DV='math90',
                     IV='Anti90', IV_range='tumble',
                     MOD='Hyper90', MOD_levels='quantiles', 
                     quantiles_IV=c(.1, .9), quantiles_MOD=c(.25, .5, .75),
                     COVARS=c('age90month','female','grade90','minority'),
                     center = FALSE, 
                     plot_type = 'interaction')	
}
}
