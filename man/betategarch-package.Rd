\name{betategarch-package}
\alias{betategarch-package}
\alias{betategarch}
\docType{package}
\title{Estimation and simulation of the 1st. order Beta-t-EGARCH model}

\description{
Estimates and simulates the 1st. order Beta-t-EGARCH model with leverage and student t conditional errors, see Harvey and Chakravarty (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge.\cr

Let y[t] denote a financial return at time t equal to

y[t] = sigma[t]*epsilon[t]

where sigma[t] is the scale or volatility, and where epsilon[t] is IID and t-distributed with df degrees of freedom. Then the 1st. order log-volatility specifiction of the Beta-t-EGARCH model can be parametrised as

lambda[t] = delta + phi1*lambda[t-1] + kappa1*u[t-1] + kappa1star*sign[-y]*(u[t-1]+1).

phi1 is the GARCH parameter (stability requires |phi1| < 1), kappa1 is the ARCH parameter, kappa1star is the leverage or volatility-asymmetry parameter, and u[t] is the score or first derivative of the log-likelihood with respect to lambda. The score u[t] is zero-mean, IID and Beta distributed, hence the name Beta-t-EGARCH model.

}

\details{
\tabular{ll}{
Package: \tab betategarch\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2012-06-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The two main functions of the package are tegarch.sim and tegarch.est. The former simulates a 1st. order Beta-t-EGARCH model whereas the latter estimates one. tegarch.recursion() and tegarch.logl() are auxiliary functions, whereas tegarch.fit() is a wrapper to tegarch.recursion. tegarch.fit() provides fitted volatility, residuals, etc., of an estimated model.
}
