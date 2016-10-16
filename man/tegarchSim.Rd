\name{tegarchSim}
\alias{tegarchSim}

\title{ Simulate from a first order Beta-Skew-t-EGARCH model }
\description{
Simulate the y series (typically interpreted as a financial return or the error in a regression) from a first order Beta-Skew-t-EGARCH model. Optionally, the conditional scale (sigma), log-scale (lambda), conditional standard deviation (stdev), dynamic components (lambdadagger in the 1-component specification, lambda1dagger and lambda2dagger in the 2-component specification), score (u) and centred innovations (epsilon) are also returned.
}
\usage{
tegarchSim(n, omega = 0, phi1 = 0.95, phi2 = 0, kappa1 = 0.01, kappa2 = 0,
  kappastar = 0, df = 10, skew = 1, lambda.initial = NULL, verbose = FALSE)
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
  \item{lambda.initial}{NULL (default) or initial value(s) of the recursion for lambda or log-volatility. If NULL then the values are chosen automatically}
  \item{verbose}{logical, TRUE or FALSE (default). If TRUE then a matrix with n rows containing y, sigma, lambda, lambdadagger, u and epsilon is returned. If FALSE then only y is returned}
}
\details{Empty}
\value{
A \code{\link{zoo}} vector of length n or a \code{\link{zoo}} matrix with n rows, depending on the value of verbose.
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2014),	'EGARCH models with fat tails, skewness and leverage'. Computational Statistics and Data Analysis 76, pp. 320-338.\cr

Sucarrat (2013), 'betategarch: Simulation, Estimation and Forecasting of First-Order Beta-Skew-t-EGARCH models'. The R Journal (Volume 5/2), pp. 137-147.\cr
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{tegarch}}, \code{\link{zoo}} }
\examples{
##1-component specification: simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05,
  df=10, skew=0.8)

##simulate the same series, but with more output (volatility, log-volatility or
##lambda, lambdadagger, u and epsilon)
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05, df=10, skew=0.8,
  verbose=TRUE)
  
##plot the simulated values:
plot(y)

##2-component specification: simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.95, phi2=0.9, kappa1=0.01, kappa2=0.05,
  kappastar=0.03, df=10, skew=0.8)
}
\keyword{Statistical Models}
