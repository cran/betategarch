\name{tegarch.fit}
\alias{tegarch.fit}

\title{
Compute the fitted values
}
\description{
Returns the fitted values of sigma, lambda, u and epsilon. Essentially a wrapper to the tegarch.recursion function.
}
\usage{
tegarch.fit(y, par, lambda.initial = NULL, c.code = TRUE)
}
\arguments{
  \item{y}{
numeric vector, typically a financial return series.
}
  \item{par}{
numeric five-dimensional vector with the estimates of delta, phi1, kappa1, kappa1star, df.
}
  \item{lambda.initial}{
initial value of the recursion for lambda or log-volatility. Default: delta/(1-phi1)
}
  \item{c.code}{
logical, TRUE (default) or FALSE. TRUE is faster since it makes use of compiled C-code
}
}
\details{Empty}

\value{
A five-column matrix with y and the fitted values of sigma, lambda, u and epsilon.
}
\references{Harvey, A. and Chakravarty, T (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge.}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{
tegarch.recursion
}
\examples{
##simulate series with 500 observations:
set.seed(123)
x <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

##estimate a 1st. order Beta-t-EGARCH model:
est <- tegarch.est(x)

##compute fitted series:
fit <- tegarch.fit(x, est$par)

##plot sigma (fitted volatility):
plot(fit[,2])

}

\keyword{Statistical Models}
