\name{predict.tegarch}
\alias{predict.tegarch}

\title{ Generate volatility forecasts \emph{n}-steps ahead }
\description{
Generates volatility forecasts from a model fitted by \code{\link{tegarch}} (i.e. a Beta-Skew-t-EGARCH model)
}
\usage{
\method{predict}{tegarch}(object, n.ahead = 1, initial.values = NULL, n.sim = 10000,
  verbose = FALSE, ...)
}
\arguments{
  \item{object}{an object of class 'tegarch'.}
  \item{n.ahead}{the number of steps ahead for which prediction is required.}
  \item{initial.values}{a vector containing the initial values of lambda and lambdadagger (lambda1dagger and lambda2dagger for 2-component models). If NULL (default) then the fitted values associated with the last return-observation are used}
  \item{n.sim}{number of simulated skew t variates.}
  \item{verbose}{logical. If FALSE (default) then only the conditional standard deviations are returned. If TRUE then also the scale is returned.}
  \item{...}{additional arguments}
}
\details{The forecast formulas of exponential ARCH models are much more complicated than those of ordinary or non-exponential ARCH models. This is particularly the case when the conditional density is skewed. The forecast formula of the conditional scale of the Beta-Skew-t-EGARCH model is not available in closed form. Accordingly, some terms (expectations involving the skewed t) are estimated numerically by means of simulation.}
\value{
A \code{zoo} object. If verbose = FALSE, then the zoo object is a vector with the forecasted conditional standard deviations. If verbose = TRUE, then the zoo object is a matrix with forecasts of both the conditional scale and the conditional standard deviation
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2013), 'EGARCH models with fat tails, skewness and leverage', forthcoming in Computational Statistics and Data Analysis. Working paper version: Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{tegarch}}, \code{\link{predict}} }
\examples{
##simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05, df=10, skew=0.8)

##estimate a 1st. order Beta-t-EGARCH model and store the output in mymod:
mymod <- tegarch(y)

#plot forecasts of volatility 1-step ahead up to 10-steps ahead:
plot(predict(mymod, n.ahead=10))
}
\keyword{Statistical Models}
