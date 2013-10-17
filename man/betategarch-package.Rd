\name{betategarch-package}
\alias{betategarch-package}
\alias{betategarch}
\docType{package}
\title{Simulation, estimation and forecasting of Beta-Skew-t-EGARCH models}

\description{
This package provides facilities for the simulation, estimation and forecasting of first order Beta-Skew-t-EGARCH models with leverage and skewed student t conditional errors, see Harvey and Sucarrat (2013).\cr

Let y[t] denote a financial return at time t equal to

y[t] = sigma[t]*epsilon[t]

where sigma[t] > 0 is the scale or volatility (generally not equal to the conditional standard deviation), and where epsilon[t] is IID and t-distributed (possibly skewed) with df degrees of freedom. Then the first order log-volatility specifiction of the one-component Beta-Skew-t-EGARCH model can be parametrised as

sigma[t] = exp(lambda[t]),\cr
lambda[t] = omega + lambdadagger,\cr
lambdadagger[t] = phi1*lambdadagger[t-1] + kappa1*u[t-1] + kappastar*sign[-y]*(u[t-1]+1).

So the scale or volatility is given by sigma[t] = exp(lambda[t]). The omega is the unconditional or long-term log-volatility, phi1 is the GARCH parameter (|phi1| < 1 implies stability), kappa1 is the ARCH parameter, kappastar is the leverage or volatility-asymmetry parameter and u[t] is the conditional score or first derivative of the log-likelihood with respect to lambda. The score u[t] is zero-mean and IID, and (u[t]+1)/(df+1) is Beta distributed when there is no skew in the conditional density of epsilon[t]. The two-component specification is given by

sigma[t] = exp(lambda[t]),\cr
lambda[t] = omega + lambda1dagger + lambda2dagger,\cr
lambda1dagger[t] = phi1*lambdadagger[t-1] + kappa1*u[t-1],\cr
lambda2dagger[t] = phi2*lambdadagger[t-1] + kappa2*u[t-1] + kappastar*sign[-y]*(u[t-1]+1).

The first component, lambda1dagger, is interpreted as the long-term component, whereas the second component, lambda2dagger, is interpreted as the short-term component.
}

\details{
\tabular{ll}{
Package: \tab betategarch\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2013-10-17\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The two main functions of the package are \code{\link{tegarchSim}} and \code{\link{tegarch}}. The first simulates a Beta-Skew-t-EGARCH models whereas the second estimates one. The second object returns an object (a list) of class 'tegarch', and a collection of methods can be applied to this class: \code{\link{coef.tegarch}}, \code{\link{fitted.tegarch}}, \code{\link{logLik.tegarch}}, \code{\link{predict.tegarch}}, \code{\link{print.tegarch}}, \code{\link{residuals.tegarch}}, \code{\link{summary.tegarch}} and \code{\link{vcov.tegarch}}. In addition, the output produced by the \code{\link{tegarchSim}} function and the \code{\link{fitted.tegarch}} and \code{\link{residuals.tegarch}} methods are of the Z's ordered observations (\code{\link{zoo}}) class, which means a range of time-series methods are available for these objects.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2013), 'EGARCH models with fat tails, skewness and leverage', forthcoming in Computational Statistics and Data Analysis. Working paper version: Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr

}
\examples{
#simulate 1000 observations from model with default parameter values:
set.seed(123)
y <- tegarchSim(500)

#estimate and store as 'mymod':
mymod <- tegarch(y)

#print estimates and standard errors:
print(mymod)

#graph of fitted volatility (conditional standard deviation):
plot(fitted(mymod))

#plot forecasts of volatility 1-step ahead up to 10-steps ahead:
plot(predict(mymod, n.ahead=10))
}
