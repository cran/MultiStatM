% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Commutators.R
\name{EliminMatr}
\alias{EliminMatr}
\title{Elimination Matrix}
\usage{
EliminMatr(d, q, useSparse = FALSE)
}
\arguments{
\item{d}{dimension of a vector x}

\item{q}{power of the Kronecker product}

\item{useSparse}{TRUE or FALSE.}
}
\value{
Elimination matrix of order \eqn{\eta_{d,q} \times d^q= {d+q-1 \choose q}}.
If \code{useSparse=TRUE} an object of the class "dgCMatrix" is produced.
}
\description{
Eliminates the duplicated/q-plicated elements  in a T-vector of multivariate moments
and cumulants.
}
\examples{
x<-c(1,2,3)
y<-kronecker(kronecker(x,x),x)
## Distinct elements of y
z<-as.matrix(EliminMatr(3,3))\%*\%y
## Restore eliminated elements in z
as.vector(QplicMatr(3,3)\%*\%z)

}
\references{
Gy. Terdik, Multivariate statistical methods - going beyond the linear,
Springer 2021.  Section  1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
}
\seealso{
Other Matrices and commutators: 
\code{\link{EliminIndx}()},
\code{\link{QplicIndx}()},
\code{\link{QplicMatr}()},
\code{\link{SymIndx}()},
\code{\link{SymMatr}()},
\code{\link{UnivMomCum}()}
}
\concept{Matrices and commutators}
