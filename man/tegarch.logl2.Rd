\name{tegarch.logl2}
\alias{tegarch.logl2}

\title{
Log-likelihood function of the two-component Beta-skew-t-EGARCH model
}
\description{
Provides the value of the log-likelihood. Auxiliary function used by tegarch.est in the estimation of the Beta-skew-t-EGARCH model. Not intended for the average user.
}

\usage{
tegarch.logl2(y, pars, lower = -Inf, upper = Inf, lambda.initial = NULL,
  logl.penalty = -1e+101, c.code = TRUE, aux = NULL)
}

\arguments{
  \item{y}{
numeric vector, typically a financial return series.
}
  \item{pars}{
numeric vector, the parameter values.
}
  \item{lower}{
numeric vector, the lower bounds used during estimation.
}
  \item{upper}{
numeric vector, the upper bounds used during estimation.
}
  \item{lambda.initial}{
NULL (default) or initial value(s) of the recursion for lambda. If NULL, then the values are chosen automatically.
}
  \item{logl.penalty}{
NULL (default) or a numeric value. If NULL, then the values are chosen automatically.
}
  \item{c.code}{
logical, TRUE (default) or FALSE. TRUE is faster since it makes use of compiled C-code.
}
  \item{aux}{
NULL (default) or a list, se code.
}
}

\details{Empty}

\value{
The value of the log-likelihood at the provided parameter values for the given data.
}
\references{
Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{
tegarch.recursion, tegarch.logl, tegarch.est
}

\keyword{Statistical Models}
