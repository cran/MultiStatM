% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation.r
\name{VarianceKurt}
\alias{VarianceKurt}
\title{Asymptotic Variance of the estimated  kurtosis vector}
\usage{
VarianceKurt(cum)
}
\arguments{
\item{cum}{The theoretical/estimated cumulants up to the 8th order in vector form}
}
\value{
The matrix of theoretical/estimated variance
}
\description{
Warning: the function requires 8! computations, for d>3, the timing required maybe large.
}
\references{
Gy.Terdik, Multivariate statistical methods - going beyond the linear,
Springer 2021. Ch. 6, formula (6.26)
}
\seealso{
Other Estimation: 
\code{\link{SampleEVSK}()},
\code{\link{SampleKurt}()},
\code{\link{SampleMomCum}()},
\code{\link{SampleSkew}()},
\code{\link{VarianceSkew}()}
}
\concept{Estimation}
