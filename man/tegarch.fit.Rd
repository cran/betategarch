\name{tegarch.fit}
\alias{tegarch.fit}
\title{
Compute the fitted values
}
\description{
Computes the fitted values. Essentially a wrapper to the tegarch.recursion (one-component) and tegarch.recursion2 (two-component) functions.
}

\usage{
tegarch.fit(y, result, lambda.initial = NULL)
}

\arguments{
  \item{y}{numeric vector, typically a financial return series.}
  \item{result}{
a list with estimation result, for example as a result of result <- tegarch.est(y)
}
  \item{lambda.initial}{
NULL (default) or initial value(s) of the recursion for lambda. If NULL, then the values are chosen automatically.
}
}

\details{Empty}

\value{
A matrix with y and the fitted values of, amongst other, sigma, the log-scale (lambda), the conditional standard deviation (stdev), u, epsilon and the standardised residuals (residstd).
}
\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Chakravarty (2008), 'Beta-t-(E)GARCH', Cambridge Working Papers in Economics 0840, Faculty of Economics, University of Cambridge. Nelson, Daniel B. (1991): 'Conditional Heteroskedasticity in Asset Returns: A New Approach', Econometrica 59, pp. 347-370.\cr

Harvey and Sucarrat (2012), ' EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{
tegarch.recursion, tegarch.recursion2
}
\examples{
##simulate series with 500 observations:
set.seed(123)
y <- tegarch.sim(500, omega=0.01, phi1=0.9, kappa1=0.1, kappastar=0.05, df=10, skew=0.8)

##estimate a 1st. order Beta-skew-t-EGARCH model:
mymod <- tegarch.est(y)

##compute fitted values:
fitvals <- tegarch.fit(y, mymod)

##make plot of the conditional standard deviaiton:
plot(fitvals[,"stdev"])

##make histogram of the standardised residuals:
hist(fitvals[,"residstd"])
}

\keyword{Statistical Models}
