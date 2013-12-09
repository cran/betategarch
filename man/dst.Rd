\name{dST}
\alias{dST}
\alias{rST}
\alias{STmean}
\alias{STvar}
\alias{STskewness}
\alias{STkurtosis}

\title{The skewed t distribution}
\description{
Density, random number generation, mean, variance, skewness and kurtosis functions for the uncentred skewed t distribution. The skewing method is that of Fernandez and Steel (1998).
}
\usage{
dST(y, df = 10, sd = 1, skew = 1, log = FALSE)
rST(n, df = 10, skew = 1)
STmean(df, skew = 1)
STvar(df, skew = 1)
STskewness(df, skew = 1)
STkurtosis(df, skew = 1)
}
\arguments{
  \item{y}{numeric vector of quantiles}
  \item{n}{integer, the number of observations}
  \item{df}{degrees of freedom, greater than 0 and less than Inf}
  \item{sd}{scale, greater than 0 and less than Inf}
  \item{skew}{skewness, greater than 0 and less than Inf. Symmetry obtains when skew = 1 (default).}
  \item{log}{logical. TRUE returns the natural log of the density value, FALSE (default) returns the density value.}
}
\details{Empty}
\value{
  \item{dST:}{a numeric value, either the density value or the natural log of the density value}
  \item{rST:}{a numeric vector with n random numbers}
  \item{STmean:}{The mean of an uncentred skewed t variable}
  \item{STvar:}{The variance of an uncentred skewed t variable}
  \item{STskewness:}{3rd. moment of a standardised skewed t variable}
  \item{STkurtosis:}{4th. moment of a standardised skewed t variable}
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr
}
\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{ \code{\link{tegarchSim}} }

\examples{
##generate 1000 random numbers from the skewed t:
set.seed(123)
eps <- rST(500, df=5) #symmetric t
eps <- rST(500, df=5, skew=0.8) #skewed to the left
eps <- rST(500, df=5, skew=2) #skewed to the right

##compare empirical mean with analytical:
mean(eps)
STmean(5, skew=2)

##compare empirical variance with analytical:
var(eps)
STvar(5, skew=2)
}
\keyword{Statistical Models}
