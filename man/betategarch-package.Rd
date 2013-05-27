\name{betategarch-package}
\alias{betategarch-package}
\alias{betategarch}
\docType{package}
\title{Estimation and simulation of Beta-t-EGARCH models}

\description{
This package provides facilities for the estimation and simulation of first order Beta-t-EGARCH models with leverage and student t conditional errors, possibly skewed, see Harvey and Sucarrat (2012), and Harvey and Chakravarty (2008).\cr

Let y[t] denote a financial return at time t equal to

y[t] = sigma[t]*epsilon[t]

where sigma[t] > 0 is the scale or volatility (generally not equal to the conditional standard deviation), and where epsilon[t] is IID and t-distributed (possibly skewed) with df degrees of freedom. Then the first order log-volatility specifiction of the one-component Beta-skew-t-EGARCH model can be parametrised as

lambda[t] = omega + lambdadagger,\cr
lambdadagger[t] = phi1*lambdadagger[t-1] + kappa1*u[t-1] + kappastar*sign[-y]*(u[t-1]+1).

The scale or volatility is given by sigma[t] = exp(lambda[t]). omega is the unconditional or long-term log-volatility, phi1 is the GARCH parameter (stability requires |phi1| < 1), kappa1 is the ARCH parameter, kappastar is the leverage or volatility-asymmetry parameter, and u[t] is the conditional score or first derivative of the log-likelihood with respect to lambda. The score u[t] is zero-mean and IID, and - when there is no skew in the conditional density - (u[t]+1)/(df+1) is Beta distributed.
}

\details{
\tabular{ll}{
Package: \tab betategarch\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-05-27\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The four main functions of the package are tegarch.sim, tegarch.sim2, tegarch.est and tegarch.fit. The first two simulate, respectively, a one-component and a two-component first order Beta-skew-t-EGARCH model. tegarch.est estimates a Beta-skew-t-EGARCH model, either with one or two components, whereas tegarch.fit generates the fitted values (i.e. volatility, residuals, etc.) of an estimated model. The rest of the functions are auxiliary functions called by the four main functions.
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Chakravarty (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge.\cr

Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr

}
