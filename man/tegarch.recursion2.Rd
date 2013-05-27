\name{tegarch.recursion2}
\alias{tegarch.recursion2}

\title{
Recursion of the two-component Beta-skew-t-EGARCH model
}
\description{
Auxiliary function used by tegarch.logl and tegarch.fit to compute the lambda series. Not intended for the average user.
}

\usage{
tegarch.recursion2(y, omega = 0.1, phi1 = 0.4, phi2 = 0.2, kappa1 = 0.05,
  kappa2 = 0.1, kappastar = 0.02, df = 10, skew = 0.6, lambda.initial = NULL,
  c.code = TRUE, verbose = FALSE, aux = NULL)
}

\arguments{
  \item{y}{numeric vector, typically a financial return series or the error of a regression.}
  \item{omega}{numeric, the value of omega}
  \item{phi1}{numeric, the value of phi1}
  \item{phi2}{numeric, the value of phi2}
  \item{kappa1}{numeric, the value of kappa1}
  \item{kappa2}{numeric, the value of kappa2}
  \item{kappastar}{numeric, the value of kappastar}
  \item{df}{numeric, the value of df (degrees of freedom)}
  \item{skew}{numeric, the value of skew (skewness parameter)}
  \item{lambda.initial}{
NULL (default) or initial value(s) of the recursion for lambda. If NULL, then the values are chosen automatically.
}
  \item{c.code}{
logical, TRUE (default) or FALSE. TRUE is faster since it makes use of compiled C-code
}
  \item{verbose}{
logical, TRUE or FALSE (default). If FALSE, then only lambda is returned. If TRUE, then a matrix with y and the fitted values of, amongst other, sigma, log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd) are returned.
}
  \item{aux}{
NULL (default) or a list, se code.
}
}

\details{Empty}

\value{
A numeric vector containing the lambda values if verbose=FALSE (default). If verbose=TRUE, then a matrix then a matrix with y and the fitted values of sigma, the log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd) are returned.
}
\references{
Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{
tegarch.fit, tegarch.logl2, tegarch.recursion
}

\keyword{Statistical Models}
