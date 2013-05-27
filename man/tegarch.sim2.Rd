\name{tegarch.sim2}
\alias{tegarch.sim2}

\title{
Simulate from a first order two-component Beta-skew-t-EGARCH model
}
\description{
Simulate the y series (typically understood as a financial return or the error in a regression) from a 1st. order 2 component Beta-skew-t-EGARCH model. Optionally, the sigma, lambda, lambdadagger, u and epsilon are also returned.
}

\usage{
tegarch.sim2(n, omega = 0, phi1 = 0.95, phi2 = 0.9, kappa1 = 0.001,
  kappa2 = 0.01, kappastar = 0.005, df = 10, skew = 0.98,
  lambda.initial = NULL, verbose = FALSE)
}

\arguments{
  \item{n}{integer, length of y (i.e. no of observations)}
  \item{omega}{numeric, the value of omega}
  \item{phi1}{numeric, the value of phi1}
  \item{phi2}{numeric, the value of phi2}
  \item{kappa1}{numeric, the value of kappa1}
  \item{kappa2}{numeric, the value of kappa2}
  \item{kappastar}{numeric, the value of kappastar}
  \item{df}{numeric, the value of df (degrees of freedom)}
  \item{skew}{numeric, the value of skew (skewness parameter}
  \item{lambda.initial}{
NULL (default) or initial value(s) of the recursion for lambda or log-volatility. If NULL then the values are chosen automatically.
}
  \item{verbose}{
logical, TRUE or FALSE (default). If TRUE then a matrix with n rows containing y, sigma, lambda, lambdadagger, u and epsilon is returned. If FALSE then only y is returned.
}
}

\details{Empty}

\value{
A numeric vector of length n or a matrix with n rows, depending on the value of verbose.
}
\references{
Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{tegarch.est}

\examples{
##simulate series with 500 observations:
set.seed(123)
y <- tegarch.sim2(500, omega=0.01, phi1=0.95, phi2=0.9, kappa1=0.01, kappa2=0.05,
  kappastar=0.03, df=10, skew=0.8)

##simulate the same series, but with more output (volatility, log-volatility or
##lambda, lambda1dagger, lambda2dagger, u and epsilon)
set.seed(123)
y <- tegarch.sim2(500, omega=0.01, phi1=0.95, phi2=0.9, kappa1=0.01, kappa2=0.05,
  kappastar=0.03, df=10, skew=0.8, verbose=TRUE)

}

\keyword{Statistical Models}
