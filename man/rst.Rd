\name{rst}
\alias{rst}

\title{ Random number generator for the skewed t distribution }

\description{ Random number generation from a skewed t distribution based on the Fernandez and Steel (1998) skewing method.}

\usage{
rst(n, df = 10, sd = 1, skew = 1)
}

\arguments{
  \item{n}{integer, the number of observations}
  \item{df}{degrees of freedom, greater than 0 and less than Inf}
  \item{sd}{scale, greater than 0 and less than Inf}
  \item{skew}{skewness, greater than 0 and less than Inf. Symmetry obtains when skew = 1 (default).}
}

\details{Empty}

\value{a numeric vector with n random numbers}

\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{tegarch.sim, tegarch.sim2, dst}

\examples{
##generate 1000 random numbers from the skewed t:
set.seed(123)
eps <- rst(1000, df=5) #symmetric t
eps <- rst(1000, df=5, skew=0.8) #skewed to the left
eps <- rst(1000, df=5, skew=2) #skewed to the right

}

\keyword{Statistical Models}
