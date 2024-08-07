% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation.r
\name{SampleGC}
\alias{SampleGC}
\title{Gram-Charlier approximation to a multivariate density}
\usage{
SampleGC(X, k = 4, cum = NULL)
}
\arguments{
\item{X}{A matrix of d-variate data}

\item{k}{the order of the approximation, by default set to 4;
(k must not be smaller than 3 or greater than 8)}

\item{cum}{if NULL (default) the cumulant vector is estimated from X.
If \code{cum} is provided no estimation of cumulants is performed.}
}
\value{
The vector of the Gram-Charlier density evaluated at X
}
\description{
Provides the truncated Gram-Charlier approximation to a multivariate density. Approximation can
be up to the first k=8 cumulants.
}
\examples{
# Gram-Charlier density approximation (k=4) of data generated from
# a bivariate skew-gaussian distribution
n<-50
alpha<-c(10,0)
omega<-diag(2)
X<-rSkewNorm(n,omega,alpha)
EC<-SampleEVSK(X)
fy4<-SampleGC(X[1:5,],cum=EC)
}
\references{
Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
Springer 2021. Section 4.7.
}
