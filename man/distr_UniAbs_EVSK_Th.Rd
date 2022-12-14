% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distributions.R
\name{distr_UniAbs_EVSK_Th}
\alias{distr_UniAbs_EVSK_Th}
\title{Moments of the modulus of the Uniform distribution on the sphere}
\usage{
distr_UniAbs_EVSK_Th(d, nCum = FALSE)
}
\arguments{
\item{d}{vector-dimension}

\item{nCum}{if it is TRUE then cumulants, skewness an kurtosis are calculated}
}
\value{
The list of the first four moments in vector form
}
\description{
Moments (up to the 4th order) of the modulus of the d-variate Uniform distribution on
the sphere on (d-1)
}
\references{
Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
Springer 2021, Lemma 5.12 p.298
}
\seealso{
Other Theoretical Moments and Cumulants: 
\code{\link{distr_CFUSN_MomCum_Th}()},
\code{\link{distr_SkewNorm_EVSK_Th}()},
\code{\link{distr_SkewNorm_MomCum_Th}()},
\code{\link{distr_Uni_EVSK_Th}()},
\code{\link{distr_Uni_MomCum_Th}()},
\code{\link{distr_ZabsM_MomCum_Th}()},
\code{\link{distr_Zabs_MomCum_Th}()}

Other Multivariate distributions: 
\code{\link{distr_CFUSN_MomCum_Th}()},
\code{\link{distr_CFUSN_Rand}()},
\code{\link{distr_CFUSSD_Rand}()},
\code{\link{distr_SkewNorm_EVSK_Th}()},
\code{\link{distr_SkewNorm_MomCum_Th}()},
\code{\link{distr_SkewNorm_Rand}()},
\code{\link{distr_Uni_EVSK_Th}()},
\code{\link{distr_Uni_MomCum_Th}()},
\code{\link{distr_Uni_Rand}()},
\code{\link{distr_ZabsM_MomCum_Th}()},
\code{\link{distr_Zabs_MomCum_Th}()}
}
\concept{Multivariate distributions}
\concept{Theoretical Moments and Cumulants}
