# 1. GramCharlier
# 2. Edgeworth
# 3. CornishFisher ## working for GC right now
# 4. supporting functions for the above (not in  index)

#' Gram-Charlier approximation to a multivariate density
#'
#' Provides the truncated Gram-Charlier approximation to a multivariate density. Approximation can
#' be up to the first k=8 cumulants according to the formula
#' \deqn{f_{\mathbf{X}}\left( \mathbf{x}\right) =\left(
#'  1+\sum_{k=3}^{\infty }\frac{1}{k!}\mathbf{B}_{k}^{\intercal }\left( 0,0,
#'  \boldsymbol{\kappa }_{\mathbf{Y},3}^{\otimes },\ldots \boldsymbol{\kappa }_{
#'  \mathbf{Y},k}^{\otimes }\right) \mathbf{H}_{k}\left( \mathbf{y}|\mathbf{I}
#'  \right) \right) \varphi \left( \mathbf{x}|\boldsymbol{\mu},\boldsymbol{\Sigma }_{
#' \mathbf{X}}\right),  }
#' where the Hermite polynomial \eqn{\mathbf{H}_{k}\left( \mathbf{y}| \mathbf{I}\right)
#' =\mathbf{H}_{k}\left( \boldsymbol{\Sigma }^{-1/2}\mathbf{x} \right) } corresponds to the \emph{standard} Gaussian variate,
#' the cumulants are the cumulants of the standardized variate \eqn{\mathbf{Y}=
#' \boldsymbol{\Sigma }^{-1/2}\left( \mathbf{X}-\boldsymbol{\mu }\right) } of
#' \eqn{\mathbf{X}}, (\eqn{\boldsymbol{\mu }={E}\mathbf{X}}) and \eqn{\varphi \left(
#'  \mathbf{x}\boldsymbol{|\boldsymbol{\mu},\Sigma }\right)} denotes the multivariate normal
#' density function with mean \eqn{\boldsymbol{\mu}} and variance matrix \eqn{\boldsymbol{
#'  \Sigma }}.
#'
#'
#' @param X A matrix of d-variate data
#' @param cum  a list containing the raw (unstandardized) cumulant vectors of X. At least the
#' first 3 cumulants need to be provided
#' @return The vector of the Gram-Charlier density evaluated at X
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021. Section 4.7.
#'
#' @examples
#' # Gram-Charlier density approximation (k=4) of data generated from
#' # a bivariate skew-gaussian distribution
#' n<-500
#' alpha<-c(10,0)
#' omega<-diag(2)
#' X<-rSkewNorm(n,omega,alpha)
#' EC<-SampleMomCum(X,r=4,centering=FALSE,scaling=FALSE)
#' EC<-EC$estCum.r  ## (estimated) raw cumulants of X
#' fx4<-GramCharlier(X[1:50,],cum=EC)
#' @family Approximations
#'
#' @export
GramCharlier<-function(X,cum){
  k=length(cum)
  if (k<3) stop("k must be greater than 2")
  if (k>8) stop("k cannot be greater than 8")
  if (is.vector(X)) stop(" X must be a data matrix")

  d<-dim(X)[[2]]
  #  cum are used to standardize
  EC<-cum
  ECz <- vector(mode="list", length=k)
  ECz[[1]]=rep(0,d)
  ECz[[2]]=diag(d)
  z1<-t(apply(X,1, function(x) x-as.vector(EC[[1]])))
  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx <- svd(cx)
  sm12 <- svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  if (is.vector(z1)) {as.matrix(z1);Z<-t(sm12%*%z1) }  else {Z<-t(sm12%*%t(z1))}
  for (j in c(3:k)){
    ECz[[j]]<-.kron_power(sm12,j)%*%EC[[j]]
  }
  phi<-mvtnorm::dmvnorm(X,EC[[1]],cx)

  HeN_L <- .SampleHermiteN_d(Z,k)

  gy<-1
  for (j in 3:min(k,5)){
    HeN <- HeN_L[[j]] #
    gy<-gy+HeN%*%ECz[[j]]/factorial(j)
  }
  if (k>5){
    if (k==6) {B6<-SymIndx(ECz[[6]]+10*kronecker(ECz[[3]],ECz[[3]]),d,6)
    Bell<-list(B6)
    }
    if (k==7) {B6<-SymIndx(ECz[[6]]+10*kronecker(ECz[[3]],ECz[[3]]),d,6)
    B7<-SymIndx(ECz[[7]]+35*kronecker(ECz[[3]],ECz[[4]]),d,7)
    Bell<-list(B6,B7)
    }
    if (k==8) {B6<-SymIndx(ECz[[6]]+10*kronecker(ECz[[3]],ECz[[3]]),d,6)
    B7<-SymIndx(ECz[[7]]+35*kronecker(ECz[[3]],ECz[[4]]),d,7)
    B8<-SymIndx(ECz[[8]]+56*kronecker(ECz[[3]],ECz[[5]])+35*kronecker(ECz[[4]],ECz[[4]]),d,8)
    Bell<-list(B6,B7,B8)
    }
    for (j in 6:min(k,8)){
      HNe <- HeN_L[[j]] #
      gy<-gy+HNe%*%Bell[[j-5]]/factorial(j)
    }
  }

  GC<-as.vector(c(gy)*phi)
  return(GC)
}




######## EDG #######################
#' Edgeworth expansion of a multivariate density
#'
#' Provides the truncated Edgeworth approximation to a multivariate density of
#' \eqn{W = \sqrt{n} \bar{X}}
#' Approximation can use up to the first k=8 cumulants. The function implements the formula
#' \deqn{f_{\mathbf{W}^{\left( n\right) }}\left( \mathbf{w}\right)  =\left(
#' 1+\sum_{k=1}^{\infty }\frac{n^{-k/2}}{k!}\mathbf{B}_{k}\left( \frac{%
#'    \boldsymbol{\kappa }_{\mathbf{Y},3}^{\otimes \intercal }\mathbf{H}_{3}\left(
#'      \mathbf{z}|\mathbf{I}\right) }{6},\ldots ,\frac{
#'        \boldsymbol{\kappa }_{\mathbf{Y},k+2}^{\otimes }\mathbf{H}_{k+2}\left(
#'          \mathbf{z}|\mathbf{I}\right) }{\left( k+1\right) \left( k+2\right) }\right)
#'  \right) \varphi \left( \mathbf{w}|\boldsymbol{\Sigma }_{\mathbf{X}}\right)}
#' where \eqn{\mathbf{z}=\boldsymbol{\Sigma }_{\mathbf{X}}^{-1/2}(\mathbf{X}-\boldsymbol{\mu}_{\mathbf{X}})},
#' \eqn{\mathbf{B}_{k}} denote the T-Bell Polynomials and \eqn{\varphi} denotes the
#' multivariate normal density. The case \eqn{n=1} provides and approximation
#' to the density of \eqn{\mathbf{X}} and can be compared to the \code{GramCharlier} approximation.
#'
#' @param X A matrix of d-variate data
#' @param cum  a list containing the raw (unstandardized) cumulant vectors of X. At least the
#' first 3 cumulants need to be provided.
#' @param n the number of terms in the mean \eqn{\bar{\mathbf{X}}}
#' @return The vector of the Edgeworth density evaluated at X
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021. Section 4.7.
#'
#' @examples
#' # Edgeworth density approximation (k=4) of data generated from
#' # a bivariate skew-gaussian distribution
#' n<-500
#' alpha<-c(10,0)
#' omega<-diag(2)
#' X<-rSkewNorm(n,omega,alpha)
#' EC<-SampleMomCum(X,r=4,centering=FALSE,scaling=FALSE)
#' EC<-EC$estCum.r  ## (estimated) raw cumulants of X
#' fx4<-Edgeworth(X[1:50,],cum=EC,n=1)
#'
#' @family Approximations
#' @export
Edgeworth<-function(X,cum,n=1) {
  k=length(cum)
  if (k<3) stop("k must be greater than 2")
  if (k>8) stop("k cannot be greater than 8")
  if (is.vector(X)) stop(" X must be a data matrix")

  d<-dim(X)[[2]]

  EC<-cum
  ECz <- vector(mode="list", length=k)
  ECz[[1]]=rep(0,d)
  ECz[[2]]=diag(d)
  z1<-t(apply(X,1, function(x) x-as.vector(EC[[1]])))
  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx <- svd(cx)
  sm12 <- svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  if (is.vector(z1)) {as.matrix(z1);Z<-t(sm12%*%z1) }  else {Z<-t(sm12%*%t(z1))}
  for (j in c(3:k)) {
    ECz[[j]]<-.kron_power(sm12,j)%*%EC[[j]]
  }
  phi<-mvtnorm::dmvnorm(X,EC[[1]],cx)

  for (j in 1:(k-2)){ po <- k+(j-1)*2}

  HeN_L <- .SampleHermiteN_d(Z,po)

  gy <- vector(mode="list", length=(k-2))

  for (j in 3:min(k,8)) {
    if (j==3) {gy[[1]] <-  (n^((j-2)/2)*factorial(j-2))^(-1) * HeN_L[[3]]%*%ECz[[3]]/6  }
    if (j==4) {gy[[2]] <- (n^((j-2)/2)*factorial(j-2))^(-1) * (HeN_L[[4]]%*%ECz[[4]]/12 + HeN_L[[6]]%*%kronecker(ECz[[3]],ECz[[3]])/36)}
    if (j==5) {
      gy[[3]] <-  (n^((j-2)/2)*factorial(j-2))^(-1)  * (
        HeN_L[[5]] %*% ECz[[5]] / 20 +
          3 * HeN_L[[7]] %*% kronecker(ECz[[3]], ECz[[4]]) / 72 +
          HeN_L[[9]] %*% kronecker(kronecker(ECz[[3]], ECz[[3]]), ECz[[3]]) / 216
      )
    }

    if (j==6) {
      gy[[4]] <- (n^((j-2)/2)*factorial(j-2))^(-1) * (
        HeN_L[[6]] %*% ECz[[6]] / 30 +
          (3 * HeN_L[[8]] %*% kronecker(ECz[[4]], ECz[[4]]) / 144 +
             4 * HeN_L[[8]] %*% kronecker(ECz[[3]], ECz[[5]]) / 120) +
          6 * HeN_L[[10]] %*% kronecker(kronecker(ECz[[3]], ECz[[3]]), ECz[[4]]) / 432 +
          HeN_L[[12]] %*% .kron_power(ECz[[3]],4) / 1296
      )
    }
    if (j==7) {
      gy[[5]] <- (1/n^((j-2)/2)*factorial(j-2)) * (
        HeN_L[[7]] %*% ECz[[7]] / 42 +
          10 * HeN_L[[9]] %*% kronecker(ECz[[4]], ECz[[5]]) / 240 +
          5 * HeN_L[[9]] %*% kronecker(ECz[[3]], ECz[[6]]) / 120 +
          10 * HeN_L[[11]] %*% kronecker(kronecker(ECz[[3]], ECz[[3]]), ECz[[5]]) / 432 +
          15 * HeN_L[[11]] %*% kronecker(kronecker(ECz[[3]], ECz[[4]]), ECz[[4]]) / 2592 +
          10 * HeN_L[[13]] %*% kronecker( .kron_power(ECz[[3]],3), ECz[[4]] ) / 864 +
          HeN_L[[15]] %*% .kron_power(ECz[[3]],5) / 7776
      )
    }

    if (j==8) {
      gy[[6]] <- (1/n^((j-2)/2)*factorial(j-2)) * (
        HeN_L[[8]] %*% ECz[[8]] / 56 +
          6 * HeN_L[[10]] %*% kronecker(ECz[[3]], ECz[[7]]) / 252 +
          15 * HeN_L[[10]] %*% kronecker(ECz[[4]], ECz[[6]]) / 360 +
          10 * HeN_L[[10]] %*% kronecker(ECz[[5]], ECz[[5]]) / 400 +
          15 * HeN_L[[12]] %*% kronecker(kronecker(ECz[[3]], ECz[[3]]), ECz[[6]]) / 1080 +
          60 * HeN_L[[12]] %*% kronecker(kronecker(ECz[[3]], ECz[[4]]), ECz[[5]]) / 1440 +
          15 * HeN_L[[12]] %*% .kron_power(ECz[[4]],3) / 1728 +
          45 * HeN_L[[14]] %*% kronecker(.kron_power(ECz[[3]],2), .kron_power(ECz[[4]],2)) / 5184 +
          20 * HeN_L[[14]] %*% kronecker(.kron_power(ECz[[3]],3), ECz[[5]]) / 4320 +
          15 * HeN_L[[16]] %*% kronecker(.kron_power(ECz[[3]],4), ECz[[4]]) / 15552 +
          HeN_L[[18]] %*% .kron_power(ECz[[3]],6) / 46656
      )
    }


  }

  fgy<- 1
  for (j in 1:(k-2)) {
    term <- gy[[j]]
    fgy <- fgy + term
  }


  Edge<-as.vector(c(fgy)*phi)
  return(Edge)
}

########################




