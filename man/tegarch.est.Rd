\name{tegarch.est}
\alias{tegarch.est}

\title{
Estimates a 1st. order Beta-t-EGARCH model
}
\description{
Provides Maximum Likelihood (ML) estimation of a 1s.t order Beta-t-EGARCH model and computes (optionally) the Hessian. EGARCH is short for Exponential Generalised Autoregressive Heteroscedasticity, a label that is often associated with Nelson's (1991) model. A shortcoming with Nelson's model is that the unconditional moments do not exist in general when the conditional density is t distributed. This is the major motivation behind Andrew Harvey's Beta-t-EGARCH model, which is an exponential GARCH model with leverage and t distributed conditional errors.
}
\usage{
tegarch.est(y, initial.values = c(0.001, 0.9, 0.02, 0.01, 10), lower = c(-Inf, -0.999999999, -Inf, -Inf, 2.000001), upper = c(+Inf, 0.999999999, Inf, Inf, Inf), lambda.initial = NULL, compute.hessian = TRUE, c.code = FALSE, na.replace = NA, verbose = TRUE, ...)
}
\arguments{
  \item{y}{
numeric vector, typically a financial return series.
}
  \item{initial.values}{
numeric five-dimensional vector with the initial values of: delta, phi1, kappa1, kappa1star, df.
}
  \item{lower}{
numeric five-dimensional vector with the lower bounds of the estimated parameter vector, see the documentation of the nlminb function.
}
  \item{upper}{
numeric five-dimensional vector with the upper bounds of the estimated parameter vector, see the documentation of the nlminb function.
}
  \item{lambda.initial}{
initial value of the recursion for lambda or log-volatility. Default: delta/(1-phi1)
}
  \item{compute.hessian}{
logical, TRUE (default) or FALSE. If TRUE then the Hessian is computed numerically via the internal optimhess function, see the documentation of the optim function.
}
  \item{c.code}{
logical, TRUE (not implemented yet) or FALSE (default). Until implemented (hopefully sometime in the near future), TRUE will return an error message.
}
  \item{na.replace}{
numeric five-dimensional vector. Sometimes estimation can result in NA values which is fatal for simulations. One way of avoiding this is to replace NAs with user-provided values (the values are sent to the log-likelihood function tegarch.logl).
}
  \item{verbose}{
logical TRUE (default) or FALSE. FALSE provides marginally faster estimation.
}
  \item{\dots}{
further arguments passed to the nlminb.
}
}

\details{
Let y[t] denote a financial return given by

y[t] = sigma[t]*epsilon[t]

where epsilon[t] is IID and t-distributed with df degrees of freedom, and where sigma[t] is a time-varying scale or volatility. The 1st. order log-volatility specifiction of the Beta-t-EGARCH model of Harvey (2011), is given by

lambda[t] = delta + phi1*lambda[t-1] + kappa1*u[t-1] + kappa1star*sign[-y]*(u[t-1]+1).

phi1 is the GARCH parameter (stability requires |phi1| < 1), kappa1 is the ARCH parameter, kappa1star is a leverage or volatility-asymmetry parameter, u[t] is the score or first derivative of the log-likelihood with respect to lambda. The score u[t] is IID and Beta distributed, hence the name Beta-t-EGARCH model.

}

\value{
Returns a list with the following elements:
  \item{par}{estimates of delta, phi1, kappa1, kappa1star and df}
  \item{objective}{the value of the log-likelihood at the maximum}
  \item{convergence}{An integer code. 0 indicates successful convergence, see the documentation of nlminb}
  \item{iterations}{Number of iterations performed, see the documentation of nlminb}
  \item{evaluations}{Number of objective function and gradient function evaluations, see the documentation of nlminb}
  \item{message}{A character string giving any additional information returned by the optimizer, or NULL. For details, see PORT documentation and the documentation of nlminb}
}
\references{
Harvey, Andrew (2011): 'The Beta-t-EGARCH Model', unpublished working paper
Nelson, Daniel B. (1991): 'Conditional Heteroskedasticity in Asset Returns: A New Approach', Econometrica 59, pp. 347-370
}
\author{
Genaro Sucarrat (http://www.sucarrat.net/)
}
\note{
}

\seealso{
tegarch.recursion, tegarch.logl
}
\examples{
##simulate series with 500 observations:
set.seed(123)
x <- tegarch.sim(500, delta=0.01, phi1=0.9, kappa1=0.1, kappa1star=0.05, df=10)

##estimate a 1st. order Beta-t-EGARCH model:
tegarch.est(x)

##estimate the same model but faster:
tegarch.est(x, compute.hessian=FALSE) #Do not compute Hessian; speeds up estimation

##estimate the same model but with non-default lower bounds (see nlminb documentation):
tegarch.est(x, lower=c(-1, 0.7, -1, -1, 4.1))

}
\keyword{}
