\name{apple}
\alias{apple}
\docType{data}
\title{Daily Apple stock returns
}
\description{The dataset contains two variables, day and appleret. Day is the date of the return and appleret is the daily (closing value) log-return in percent of the Apple stock over the period 10 September 1985 - 10 May 2011 (a total of 6835 observations).}
\usage{data(apple)}
\format{
  A data frame with 6835 observations:
  \describe{
    \item{\code{day}}{a factor}
    \item{\code{appleret}}{a numeric vector}
  }
}
\details{The data is studied in more detail in Harvey and Sucarrat (2012).}
\source{The source of the original raw data is http://yahoo.finance.com/.}
\references{
Harvey and Sucarrat (2012), 'EGARCH models with fat tails, skewness and leverage', Cambridge Working Papers in Economics 1236, Faculty of Economics, University of Cambridge.\cr
}
\examples{
data(apple) #load data into workspace
tegarch.est(apple[,"appleret"]) #estimate volatility model of Apple returns
}
\keyword{datasets}
