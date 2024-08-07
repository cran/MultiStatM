% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CallFcts.R
\name{SampleSkew}
\alias{SampleSkew}
\title{Estimation of Sample Skewness (Mardia, MRSz)}
\usage{
SampleSkew(x, Type = c("Mardia", "MRSz"))
}
\arguments{
\item{x}{A matrix of multivariate data.}

\item{Type}{A character string specifying the type of skewness index to estimate. Use "Mardia" for Mardia's skewness index or "MRSz" for the Mori-Rohatgi-Szekely skewness vector and index.}
}
\value{
A list containing the estimated skewness index or vector and the associated p-value under the Gaussian hypothesis.
\item{Mardia.Skewness}{The skewness index when \code{Type} is "Mardia".}
\item{MRSz.Skewness.Vector}{The skewness vector when \code{Type} is "MRSz".}
\item{MRSz.Skewness.Index}{The skewness index when \code{Type} is "MRSz".}
\item{p.value}{The p-value under the Gaussian hypothesis for the estimated skewness.}
}
\description{
Estimates the sample skewness index based on the specified method: Mardia or MRSz.
}
\examples{
# Mardia's skewness example
x <- matrix(rnorm(100*5), ncol=5)
SampleSkew(x, Type = "Mardia")

# MRSz's skewness example
SampleSkew(x, Type = "MRSz")

}
\references{
Gy.Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Example 6.1 and 6.2.

S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate skewness and kurtosis. Sankhya A, 83(2), 607-644.
}
\seealso{
Other Estimation: 
\code{\link{SampleEVSK}()},
\code{\link{SampleKurt}()},
\code{\link{SampleMomCum}()},
\code{\link{VarianceKurt}()},
\code{\link{VarianceSkew}()}
}
\concept{Estimation}
