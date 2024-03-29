% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Commutators.R
\name{indx_Commutator_Kperm}
\alias{indx_Commutator_Kperm}
\title{Index vector for commutation of T-products of any number of vectors}
\usage{
indx_Commutator_Kperm(perm, dims)
}
\arguments{
\item{perm}{vector indicating the permutation of the order
in the Kronecker product,}

\item{dims}{vector indicating the dimensions of the vectors,
use dims <- d if all dimensions are equal}
}
\value{
An index vector to produce the permutation
}
\description{
Produces any permutation of kronecker products of vectors of any length.
Same results as matr_Commutator_Kperm.
}
\examples{
a1<-c(1,2)
a2<-c(2,3,4)
a3<-c(1,3)
p1<-a1\%x\%a2\%x\%a3
p1[indx_Commutator_Kperm(c(3,1,2),c(2,3,2))]
## Same as
a3\%x\%a1\%x\%a2
## Same as
as.vector(matr_Commutator_Kperm(c(3,1,2),c(2,3,2))\%*\%p1)
}
\references{
Holmquist B (1996) The d-variate vector Hermite polynomial of order. Linear Algebra
and its Applications 237/238, 155-190.

Gy., Terdik, Multivariate statistical methods - going beyond the linear,
Springer 2021, 1.2.4 Commuting T-Products of Vectors.
}
\seealso{
Other Matrices and commutators: 
\code{\link{indx_Commutator_Kmn}()},
\code{\link{indx_Commutator_Mixing}()},
\code{\link{indx_Commutator_Moment}()},
\code{\link{indx_Elimination}()},
\code{\link{indx_Qplication}()},
\code{\link{indx_Symmetry}()},
\code{\link{indx_UnivMomCum}()},
\code{\link{matr_Commutator_Kmn}()},
\code{\link{matr_Commutator_Kperm}()},
\code{\link{matr_Commutator_Mixing}()},
\code{\link{matr_Commutator_Moment}()},
\code{\link{matr_Elimination}()},
\code{\link{matr_Qplication}()},
\code{\link{matr_Symmetry}()}
}
\concept{Matrices and commutators}
