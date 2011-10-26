\name{tegarch.sim}
\alias{tegarch.sim}

\title{
Simulate from a 1st. order Beta-t-EGARCH model
}
\description{
Generate the y series from a 1st. order Beta-t-EGARCH model and, optionally, the sigma, lambda, u and epsilon.
}
\usage{
tegarch.sim(n, delta = 0, phi1 = 0.95, kappa1 = 0.05, kappa1star = 0, df = 10,
  lambda.initial = NULL, verbose = FALSE)
}

\arguments{
  \item{n}{
integer, length of y
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
numeric, the value of df (degrees of freedom parameter of the t distribution)
}
  \item{lambda.initial}{
numeric, the initial value of the recursion for lambda or log-volatility. Default: delta/(1-phi1)
}
  \item{verbose}{
logical, TRUE or FALSE (default). If TRUE then a five-column matrix with n rows containing y, sigma, lambd, u and epsilon is returned. If FALSE then only y is returned.
}
}
\details{
}
\value{
A numeric vector of length n or a five-column matrix with n rows, depending on the value of verbose.
}
\references{
}
\author{
Genaro Sucarrat, http://www.sucarrat.net/
}
\note{
}


\seealso{
}
\examples{
##simulate series with 500 observations:
set.seed(123)
x <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

##simulate the same series, but with more output (volatility, log-volatility or lambda, u and epsilon)
set.seed(123)
y <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10, verbose=TRUE)

}

\keyword{}
