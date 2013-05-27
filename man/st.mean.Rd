\name{st.mean}
\alias{st.mean}

\title{ Compute mean of a skewed t variable }

\description{ This function computes the mean of a skewed t variable. }

\usage{
st.mean(df = 10, sd = 1, skew = 1)
}

\arguments{
  \item{df}{degrees of freedom, greater than 0 and less than Inf}
  \item{sd}{scale, greater than 0 and less than Inf}
  \item{skew}{skewness, greater than 0 and less than Inf. Symmetry obtains when skew = 1 (default).}
}

\details{Empty}

\value{The mean of skewed t variable}

\references{
Fernandez and Steel (1998), 'On Bayesian Modeling of Fat Tails and Skewness', Journal of the American Statistical Association 93, pp. 359-371.\cr

Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}

\author{Genaro Sucarrat, http://www.sucarrat.net/}

\note{Empty}

\seealso{rst, dst}

\keyword{Statistical Models}
