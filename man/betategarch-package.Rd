\name{betategarch-package}
\alias{betategarch-package}
\alias{betategarch}
\docType{package}
\title{
Estimation and simulation of the 1st. order Beta-t-EGARCH model
}
\description{
Estimates and simulates the 1st. order Beta-t-EGARCH model of Andrew Harvey (2011) with leverage and student t conditional errors. If y[t] denotes a financial return equal to

y[t] = sigma[t]*epsilon[t]

where epsilon[t] is IID and t-distributed with df degrees of freedom, and where sigma[t] is a time-varying scale or volatility, then the 1st. order log-volatility specifiction of the Beta-t-EGARCH model of Harvey (2011) is given by

lambda[t] = delta + phi1*lambda[t-1] + kappa1*u[t-1] + kappa1star*sign[-y]*(u[t-1]+1).

phi1 is the GARCH parameter (stability requires |phi1| < 1), kappa1 is the ARCH parameter, kappa1star is the leverage or volatility-asymmetry parameter, and u[t] is the score or first derivative of the log-likelihood with respect to lambda. The score u[t] is IID and Beta distributed, hence the name Beta-t-EGARCH model.

}
\details{
\tabular{ll}{
Package: \tab betategarch\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-27\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

The two main functions of the package are tegarch.est and tegarch.sim. The former estimates a 1st. order Beta-t-EGARCH model whereas the latter simulates one. tegarch.recursion() and tegarch.logl() are auxiliary functions, whereas tegarch.fit() is a wrapper to tegarch.recursion and provides fitted volatility, residuals, etc., of an estimated model.

}
