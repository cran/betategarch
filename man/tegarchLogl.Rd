\name{tegarchLogl}
\alias{tegarchLogl}
\alias{tegarchLogl2}
\alias{tegarchRecursion}
\alias{tegarchRecursion2}

\title{ Auxiliary functions }
\description{
tegarchLogl, tegarchLogl2, tegarchRecursion and tegarchRecursion2 are auxiliary functions called by \code{\link{tegarch}}, and which are not intended to be used for the average user. Henceforth they are thusonly scarcely documented, but most should either be self-explanatory (for the non-average user!) or more or less documented in relation with the \code{\link{tegarch}} and \code{\link{tegarchSim}} functions.
}
\usage{
##the '2' relates to the 2-component specification:
tegarchLogl(y, pars, lower = -Inf, upper = Inf, lambda.initial = NULL,
  logl.penalty = -1e+100, c.code = TRUE, aux = NULL)
tegarchLogl2(y, pars, lower = -Inf, upper = Inf, lambda.initial = NULL,
  logl.penalty = -1e+101, c.code = TRUE, aux = NULL)
tegarchRecursion(y, omega = 0.1, phi1 = 0.4, kappa1 = 0.2, kappastar = 0.1,
  df = 10, skew = 0.6, lambda.initial = NULL, c.code = TRUE, verbose = FALSE,
  aux = NULL)
tegarchRecursion2(y, omega = 0.1, phi1 = 0.4, phi2 = 0.2, kappa1 = 0.05,
  kappa2 = 0.1, kappastar = 0.02, df = 10, skew = 0.6, lambda.initial = NULL,
  c.code = TRUE, verbose = FALSE, aux = NULL)
}
\arguments{
  \item{y}{numeric vector, typically a financial return series}
  \item{omega}{numeric}
  \item{phi1}{numeric, must be less than 1 in absolute value}
  \item{phi2}{numeric, must be less than 1 in absolute value}
  \item{kappa1}{numeric}
  \item{kappa2}{numeric}
  \item{kappastar}{numeric}
  \item{df}{numeric, the value of df (degrees of freedom)}
  \item{skew}{numeric (positive), the value of skew (skewness parameter)}
  \item{verbose}{logical. If FALSE (default) then only lambda is returned. If TRUE then a matrix with y and the fitted values of, amongst other, sigma, the log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd) are returned}
  \item{pars}{numeric vector, the parameter values}
  \item{lower}{numeric vector, the lower bounds used during estimation}
  \item{upper}{numeric vector, the upper bounds used during estimation}
  \item{lambda.initial}{NULL (default) or initial value(s) of the recursion for lambda. If NULL, then the values are chosen automatically}
  \item{logl.penalty}{numeric value}
  \item{c.code}{logical. TRUE (default) is faster since it makes use of compiled C-code}
  \item{aux}{NULL (default) or a list, se \code{\link{tegarch}} code}
}
\details{
tegarchLogl and tegarchLogl2 return the value of the log-likelihood for a 1-component and 2-component model, respectively. 
}
\value{
  \item{tegarchLogl:}{The log-likelihood value (i.e. a numeric) of a 1-component specification}
  \item{tegarchLogl2:}{The log-likelihood value (i.e. a numeric) of a 2-component specification}
  \item{tegarchRecursion:}{A numeric vector containing the lambda values if verbose=FALSE (default). If verbose=TRUE then a matrix then a matrix with y and the fitted values of sigma, the log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd) are returned}
  \item{tegarchRecursion2:}{A numeric vector containing the lambda values if verbose=FALSE (default). If verbose=TRUE, then a matrix then a matrix with y and the fitted values of sigma, the log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd) are returned}
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2014),	'EGARCH models with fat tails, skewness and leverage'. Computational Statistics and Data Analysis 76, pp. 320-338.\cr

Sucarrat (2013), 'betategarch: Simulation, Estimation and Forecasting of First-Order Beta-Skew-t-EGARCH models'. The R Journal (Volume 5/2), pp. 137-147.\cr
}
\author{Genaro Sucarrat, http://www.sucarrat.net/}
\seealso{
\code{\link{tegarch}}, \code{\link{tegarchSim}}, \code{\link{fitted.tegarch}} 
}
\keyword{Statistical Models}
