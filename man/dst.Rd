\name{dst}
\alias{dst}

\title{ Skewed t density}

\description{ A skewed t density based on the Fernandez and Steel (1998) skewing method. }

\usage{
dst(y, df = 10, sd = 1, skew = 1, log = FALSE)
}

\arguments{
  \item{y}{numeric vector of quantiles}
  \item{df}{degrees of freedom, greater than 0 and less than Inf}
  \item{sd}{scale, greater than 0 and less than Inf}
  \item{skew}{skewness, greater than 0 and less than Inf. Symmetry obtains when skew = 1 (default).}
  \item{log}{logical, TRUE or FALSE (default). TRUE returns the natural log of the density value, FALSE returns the density value.}
}

\details{Empty}

\value{
A numeric value, either the density value or the natural log of the density value.
}

\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{rst, st.mean}

\keyword{Statistical Models}
