\name{tegarch}
\alias{tegarch}

\title{ Estimate first order Beta-Skew-t-EGARCH models }

\description{
Fits a first order Beta-Skew-t-EGARCH model to a univariate time-series by exact Maximum Likelihood (ML) estimation. Estimation is via the \code{\link{nlminb}} function
}

\usage{
tegarch(y, asym = TRUE, skew = TRUE, components = 1, initial.values = NULL,
  lower = NULL, upper = NULL, hessian = TRUE, lambda.initial = NULL,
  c.code = TRUE, logl.penalty = NULL, aux = NULL, ...)
}

\arguments{
  \item{y}{numeric vector, typically a financial return series.}
  \item{asym}{logical. TRUE (default) includes leverage or volatility asymmetry in the log-scale specification}
  \item{skew}{logical. TRUE (default) enables and estimates the skewness in conditional density (epsilon). The skewness method is that of Fernandez and Steel (1998)}
  \item{components}{Numeric value, either 1 (default) or 2. The former estimates a 1-component model, the latter a 2-component model}
  \item{initial.values}{NULL (default) or a vector with the initial values. If NULL, then the values are automatically chosen according to model (with or without skewness, 1 or 2 components, etc.)}
  \item{lower}{NULL (default) or a vector with the lower bounds of the parameter space. If NULL, then the values are automatically chosen}
  \item{upper}{NULL (default) or a vector with the upper bounds of the parameter space. If NULL, then the values are automatically chosen}
  \item{hessian}{logical. If TRUE (default) then the Hessian is computed numerically via the optimHess function}
  \item{lambda.initial}{NULL (default) or a vector with the initial value(s) of the recursion for lambda and lambdadagger. If NULL then the values are chosen automatically}
  \item{c.code}{logical. TRUE (default) is faster since it makes use of compiled C-code}
  \item{logl.penalty}{NULL (default) or a numeric value. If NULL then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in NA and/or +/-Inf values, which are fatal for simulations. The value logl.penalty is the value returned by the log-likelihood function in the presence of NA or +/-Inf values}
  \item{aux}{NULL (default) or a list, se code. Useful for simulations (speeds them up)}
  \item{\dots}{further arguments passed to the nlminb function}
}

\value{
Returns a list of class 'tegarch' with the following elements:
  \item{y}{the series used for estimation.}
  \item{date}{date and time of estimation.}
  \item{initial.values}{initial values used in estimation.}
  \item{lower}{lower bounds used in estimation.}
  \item{upper}{upper bounds used in estimation.}
  \item{lambda.initial}{initial values of lambda provided by the user, if any.}
  \item{model}{type of model estimated.}
  \item{hessian}{the numerically estimated Hessian.}
  \item{sic}{the value of the Schwarz (1978) information criterion.}
  \item{par}{parameter estimates.}
  \item{objective}{value of the log-likelihood at the maximum.}
  \item{convergence}{an integer code. 0 indicates successful convergence, see the documentation of nlminb.}
  \item{iterations}{number of iterations, see the documentation of nlminb.}
  \item{evaluations}{number of evaluations of the objective and gradient functions, see the documentation of nlminb.}
  \item{message}{a character string giving any additional information returned by the optimizer, or NULL. For details, see PORT documentation and the nlminb documentation.}
  \item{NOTE}{an additional message returned if one tries to estimate a 2-component model without leverage.}
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Nelson, Daniel B. (1991): 'Conditional Heteroskedasticity in Asset Returns: A New Approach', Econometrica 59, pp. 347-370.\cr

Harvey and Sucarrat (2013), 'EGARCH models with fat tails, skewness and leverage', forthcoming in Computational Statistics and Data Analysis. Working paper version: Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr

Schwarz (1978), 'Estimating the Dimension of a Model', The Annals of Statistics 6,  pp. 461-464.\cr
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{tegarchSim}}, \code{\link{coef.tegarch}}, \code{\link{fitted.tegarch}}, \code{\link{logLik.tegarch}}, \code{\link{predict.tegarch}}, \code{\link{print.tegarch}}, \code{\link{residuals.tegarch}}, \code{\link{summary.tegarch}}, \code{\link{vcov.tegarch}}
}

\examples{
##simulate series with 500 observations:
set.seed(123)
y <- tegarchSim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05,
  df=10, skew=0.8)

##estimate a 1st. order Beta-t-EGARCH model and store the output in mymod:
mymod <- tegarch(y)

#print estimates and standard errors:
print(mymod)

#graph of fitted volatility (conditional standard deviation):
plot(fitted(mymod))

#graph of fitted volatility and more:
plot(fitted(mymod, verbose=TRUE))

#plot forecasts of volatility 1-step ahead up to 20-steps ahead:
plot(predict(mymod, n.ahead=20))

#full variance-covariance matrix:
vcov(mymod)
}
\keyword{Statistical Models}
