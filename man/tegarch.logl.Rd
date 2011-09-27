\name{tegarch.logl}
\alias{tegarch.logl}
\title{
Log-likelihood function
}
\description{
Provides the value of the log-likelihood. Auxiliary function used by tegarch.est in the estimation of the Beta-t-EGARCH model.
}
\usage{
tegarch.logl(y, delta = 0, phi1 = 0.8, kappa1 = 0.1, kappa1star = 0.01, df = 10, lower = c(-Inf, -0.999999999, -Inf, -Inf, 2.1), upper = c(Inf, 0.999999999, Inf, Inf, Inf), lambda.initial = NULL, c.code = FALSE, na.replace = rep(NA,5))
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
  \item{lower}{
numeric five-dimensional vector with the lower bounds of the estimated parameter vector, see the documentation of the nlminb function
}
  \item{upper}{
numeric five-dimensional vector with the upper bounds of the estimated parameter vector, see the documentation of the nlminb function
}
  \item{lambda.initial}{
initial value of the recursion for lambda or log-volatility. Default: delta/(1-phi1)
}
  \item{c.code}{
logical, TRUE (not implemented yet) or FALSE (default). Until implemented (hopefully sometime in the near future), TRUE will return an error message
}
  \item{na.replace}{
numeric five-dimensional vector. Sometimes estimation can result in NA values which is fatal for simulations. One way of avoiding this is to replace NAs with user-provided values.
}
}
\details{
}
\value{
Numeric value, the value of the log-likelihood at the provided parameter values and given the data
}
\references{
}
\author{
Genaro Sucarrat, http://www.sucarrat.net/
}
\note{
}

\seealso{
tegarch.recursion
}
\examples{
##simulate series with 500 observations:
set.seed(123)
x <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

##compute the associated log-likelihood:
tegarch.logl(x, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

}

\keyword{}
