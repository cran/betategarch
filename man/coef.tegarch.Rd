\name{coef.tegarch}
\alias{coef.tegarch}
\alias{fitted.tegarch}
\alias{logLik.tegarch}
\alias{print.tegarch}
\alias{residuals.tegarch}
\alias{summary.tegarch}
\alias{vcov.tegarch}

\title{ Extraction methods for 'tegarch' objects }
\description{
Extraction methods for objects of class 'tegarch' (i.e. the result of estimating a Beta-Skew-t-EGARCH model)
}
\usage{
\method{coef}{tegarch}(object, ...)
\method{fitted}{tegarch}(object, verbose = FALSE, ...)
\method{logLik}{tegarch}(object, ...)
\method{print}{tegarch}(x, ...)
\method{residuals}{tegarch}(object, standardised = TRUE, ...)
\method{summary}{tegarch}(object, verbose = FALSE, ...)
\method{vcov}{tegarch}(object, ...)
}
\arguments{
  \item{object}{an object of class 'tegarch'}
  \item{x}{an object of class 'tegarch'}
  \item{verbose}{logical. If FALSE (default) then only basic information is returned}
  \item{standardised}{logical. If TRUE (default) then the standardised residuals are returned. If FALSE then the scaled (by sigma) residuals are returned}
  \item{...}{additional arguments}
}
\details{Empty}
\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A zoo object. If verbose=FALSE (default), then the zoo object is a vector containing the fitted conditional standard deviations. If verbose = TRUE, then the zoo object is a matrix containing the return series y, fitted scale (sigma), fitted conditional standard deviation (stdev), fitted log-scale (lambda), dynamic component(s) (lambdadagger in the 1-component specification, lambda1dagger and lambda2dagger in the 2-compoment specification), the score (u), scaled residuals (epsilon) and standardised residuals (residstd)}
  \item{logLik:}{The value of the log-likelihood at the maximum}
  \item{print:}{Prints the most important parts of the estimation results}
  \item{residuals:}{A zoo object. If standardised = TRUE (default), then the zoo object is a vector with the standardised residuals. If standardised = FALSE, then the zoo vector contains the scaled residuals} 
  \item{summary:}{A list. If verbose = FALSE, then only the most important entries are returned. If verbose = TRUE, then all entries apart from the 1st. (the y series) is returned}
  \item{vcov:}{The variance-covariance matrix of the estimated coefficents. The matrix is obtained by inverting the numerically estimated Hessian}
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2014),	'EGARCH models with fat tails, skewness and leverage'. Computational Statistics and Data Analysis 76, pp. 320-338.\cr

Sucarrat (2013), 'betategarch: Simulation, Estimation and Forecasting of First-Order Beta-Skew-t-EGARCH models'. The R Journal (Volume 5/2), pp. 137-147.\cr

}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{tegarch}}, \code{\link{coef}}, \code{\link{fitted}}, \code{\link{logLik}}, \code{\link{predict}}, \code{\link{predict.tegarch}}, \code{\link{print}}, \code{\link{summary}}, \code{\link{vcov}} }
\examples{
#simulate 500 observations from model with default parameter values:
set.seed(123)
y <- tegarchSim(500)

#estimate and store as 'mymodel':
mymod <- tegarch(y)

#print estimation result:
print(mymod)

#extract coefficients:
coef(mymod)

#extract log-likelihood:
logLik(mymod)

#plot fitted conditional standard deviations:
plot(fitted(mymod))

#plot all the fitted series:
plot(fitted(mymod, verbose=TRUE))

#histogram of standardised residuals:
hist(residuals(mymod))
}
\keyword{Statistical Models}
