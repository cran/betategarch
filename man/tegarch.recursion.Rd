\name{tegarch.recursion}
\alias{tegarch.recursion}

\title{
Recursion of the Beta-t-EGARCH
}
\description{
To compute the log-likelihood a recursion is needed in order to generate the lambda series. This function does this
}
\usage{
tegarch.recursion(y, delta = 0.01, phi1 = 0.95, kappa1 = 0.05, kappa1star = 0,
  df = 10, lambda.initial = NULL, c.code = TRUE, verbose = FALSE)
}
\arguments{
  \item{y}{
numeric vector, typically a financial return series.
}
  \item{delta}{
numeric, the value of delta
}
  \item{phi1}{
numeric, the value of phi1
}
  \item{kappa1}{
numeric, the value of kappa1
}
  \item{kappa1star}{
numeric, the value of kappa1star
}
  \item{df}{
numeric, the value of df
}
  \item{lambda.initial}{
numeric, the initial value of the recursion for lambda or log-volatility. Default: delta/(1-phi1)
}
  \item{c.code}{
logical, TRUE (default) or FALSE. TRUE is faster since it makes use of compiled C-code
}
  \item{verbose}{
logical, TRUE or FALSE (default). If FALSE then only the lambda series in the form of a numeric vector is returned. If TRUE then a five-column matrix with y, sigma, lambda, u and epsilon is returned.
}
}
\details{Empty}

\value{
A numeric vector containing the lambda values if verbose=FALSE (default), and a five-column matrix if verbose=TRUE.
}
\references{Harvey, A. and Chakravarty, T (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge.}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{
tegarch.fit
}
\examples{
##simulate series with 500 observations:
set.seed(123)
x <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

##values of lambda (log-volatility):
xfit <- tegarch.recursion(x, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05,
  df=10)

##values of y, sigma (volatility), lambda (log-volatility), u and epsilon:
xfitmore <- tegarch.recursion(x, verbose=TRUE, delta=0.01, phi1=0.9, kappa1=0.1,
  kappa1star=0.05, df=10)

}

\keyword{Statistical Models}
