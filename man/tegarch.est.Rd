\name{tegarch.est}
\alias{tegarch.est}

\title{ Estimates a first order Beta-skew-t-EGARCH model }

\description{
Provides Maximum Likelihood (ML) estimation of a first order Beta-skew-t-EGARCH model and computes, optionally, the Hessian. Estimation is via the nlminb function. EGARCH is short for Exponential Generalised Autoregressive Heteroscedasticity, a label that is often associated with Nelson's (1991) model. A shortcoming with Nelson's model is that the unconditional moments do not exist in general when the conditional density is t distributed. This is the major motivation behind Andrew Harvey's Beta-t-EGARCH class of models, which is an exponential GARCH model with leverage and t distributed conditional errors (possibly skewed).
}

\usage{
tegarch.est(y, asym = TRUE, skew = TRUE, components = 1, initial.values = NULL,
  lower = NULL, upper = NULL, compute.hessian = FALSE, lambda.initial = NULL,
  c.code = TRUE, logl.penalty = -1e+100, verbose = TRUE, aux = NULL, ...)
}

\arguments{
  \item{y}{numeric vector, typically a financial return series.}
  \item{asym}{
logical, TRUE (default) or FALSE. TRUE includes leverage or volatility asymmetry in the log-scale specification
}
  \item{skew}{
logical, TRUE (default) or FALSE. TRUE enables and estimates the skewness in conditional density (epsilon). The skewness method is that of Fernandez and Steel (1998).
}
  \item{components}{
Numeric value, either 1 (default) or 2. The former estimates a one-component model, the latter a two-component model.
}
  \item{initial.values}{
NULL (default) or a vector with the initial values. If NULL, then then function automatically chooses the initial values according to model (with or without skewness, 1 or 2 components, etc.).
}
  \item{lower}{
NULL (default) or a vector with the lower bounds. If NULL, then the values are automatically chosen.
}
  \item{upper}{
NULL (default) or a vector with the upper bounds. If NULL, then the values are automatically chosen.
}
  \item{compute.hessian}{
logical, TRUE or FALSE (default). If TRUE then the Hessian is computed numerically using the optimHess function.
}
  \item{lambda.initial}{
NULL (default) or initial value(s) of the recursion for lambda or log-volatility. If NULL then the values are chosen automatically.
}
  \item{c.code}{
logical, TRUE (default) or FALSE. TRUE is faster since it makes use of compiled C-code
}
  \item{logl.penalty}{
NULL (default) or a numeric value. If NULL then the function chooses the value automatically. Sometimes estimation can result in NA or +/-Inf values, which can be fatal for simulations. The value logl.penalty is the value returned by the log-likelihood function in the presence of NA or +/-Inf values.
}
  \item{verbose}{
logical, TRUE (default) or FALSE. FALSE provides faster estimation, because less results-output is generated.
}
  \item{aux}{
NULL (default) or a list, se code. Useful for simulations (speeds them up).
}
  \item{\dots}{
further arguments passed to the nlminb function.
}
}

\value{
Returns a list with the following elements:
  \item{date}{date and time of estimation.}
  \item{initial.values}{initial values used in estimation.}
  \item{lower}{lower bounds used in estimation.}
  \item{upper}{upper bounds used in estimation.}
  \item{model}{type of model estimated.}
  \item{hessian.numerical}{Hessian, numerically estimated.}
  \item{par}{parameter estimates.}
  \item{objective}{value of the log-likelihood at the maximum.}
  \item{convergence}{An integer code. 0 indicates successful convergence, see the documentation of nlminb.}
  \item{iterations}{Number of iterations, see the documentation of nlminb.}
  \item{evaluations}{Number of evaluations of the objective and gradient functions, see the documentation of nlminb.}
  \item{message}{A character string giving any additional information returned by the optimizer, or NULL. For details, see PORT documentation and the nlminb documentation.}
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Chakravarty (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge.\cr

Nelson, Daniel B. (1991): 'Conditional Heteroskedasticity in Asset Returns: A New Approach', Econometrica 59, pp. 347-370.\cr

Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{tegarch.sim, tegarch.sim2, tegarch.fit}

\examples{
##simulate series with 500 observations:
set.seed(123)
y <- tegarch.sim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05, df=10, skew=0.8)

##estimate a 1st. order Beta-t-EGARCH model and store the output in mymod:
mymod <- tegarch.est(y)

##estimate the same model but with Hessian:
tegarch.est(y, compute.hessian=TRUE)

}
\keyword{Statistical Models}
