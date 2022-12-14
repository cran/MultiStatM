% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distributions.R
\name{distr_CFUSN_MomCum_Th}
\alias{distr_CFUSN_MomCum_Th}
\title{Moments and cumulants CFUSN}
\usage{
distr_CFUSN_MomCum_Th(r, d, p, Delta, nMu = FALSE)
}
\arguments{
\item{r}{The highest cumulant order}

\item{d}{The multivariate dimension and number of rows of the skewness matrix Delta}

\item{p}{The number of cols of the skewness matrix Delta}

\item{Delta}{The skewness matrix}

\item{nMu}{If set to TRUE, the list of the first r d-variate moments is provided}
}
\value{
The list of theoretical  cumulants in vector form
}
\description{
Provides the theoretical cumulants of the multivariate Canonical Fundamental
Skew Normal distribution
}
\examples{
r <- 4; d <- 2; p <- 3
Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
ieg<- eigen(diag(p)+t(Lamd)\%*\%Lamd)
V <- ieg$vectors
Delta <-Lamd \%*\% V \%*\% diag(1/sqrt(ieg$values)) \%*\% t(V)
MomCum_CFUSN <- distr_CFUSN_MomCum_Th (r,d,p,Delta)
}
\references{
Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
Springer 2021, Lemma 5.3 p.251
}
\seealso{
Other Theoretical Moments and Cumulants: 
\code{\link{distr_SkewNorm_EVSK_Th}()},
\code{\link{distr_SkewNorm_MomCum_Th}()},
\code{\link{distr_UniAbs_EVSK_Th}()},
\code{\link{distr_Uni_EVSK_Th}()},
\code{\link{distr_Uni_MomCum_Th}()},
\code{\link{distr_ZabsM_MomCum_Th}()},
\code{\link{distr_Zabs_MomCum_Th}()}

Other Multivariate distributions: 
\code{\link{distr_CFUSN_Rand}()},
\code{\link{distr_CFUSSD_Rand}()},
\code{\link{distr_SkewNorm_EVSK_Th}()},
\code{\link{distr_SkewNorm_MomCum_Th}()},
\code{\link{distr_SkewNorm_Rand}()},
\code{\link{distr_UniAbs_EVSK_Th}()},
\code{\link{distr_Uni_EVSK_Th}()},
\code{\link{distr_Uni_MomCum_Th}()},
\code{\link{distr_Uni_Rand}()},
\code{\link{distr_ZabsM_MomCum_Th}()},
\code{\link{distr_Zabs_MomCum_Th}()}
}
\concept{Multivariate distributions}
\concept{Theoretical Moments and Cumulants}
