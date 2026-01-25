###############################
# Compute Recursion, GC and Edg based on recursion
# -----
# 1. IntHermiteN - integrates up to hermite of order N
# 2 .Funct_CF     - recursive formula in (-infty,y)
# 3 .Funct_CF_low - recursive formula in (y, infty)

# 4  IntGramCharlier(x, cum , type = c("lower", "upper"))
# 5  IntEdgeworth <- function(x, cum , n = 1, type = c("lower", "upper"))

# 6 .CornishFisherM (x,cum)
# 7 .CornishFisherM_low (x,cum)

# 8 .CorFisEdgeM_low (x, cum, n)
# 9 .CorFisEdgeM (x,cum,n)

# 10  .CornishFisherM_in (y_row,cumY)
# 11  .CornishFisherM_in_low (y_row,cumY)

# 12 .CorFisEdgeM_in_low (y_row,cumY,n)
# 13 .CorFisEdgeM_in(y_row,cumY,n)

# 14 .X2Ystand (x,cum)
# 15 .Coef_GC  (cum)
# 16 .Coef_Ed (cum,n)



#############################################*
.X2Ystand <- function(x,cum){
  ########## transform to standardized x for Hermite
  ########## and its coefficients
  # x matrix n x d # x <- Xt
  # cum list   # cum <- cumt
  d <-  dim(x)[2]
  k <- length(cum)
  EC<-cum
  ECz <- vector(mode="list", length=length(cum))
  ECz[[1]]=rep(0,d)
  ECz[[2]]=diag(d)

  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx <- svd(cx)
  sm12 <- svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)

  for (j in c(3:k)){
    ECz[[j]]<-.kron_power(sm12,j)%*%EC[[j]]
  }
  z1<-t(apply(x,1, function(x) (x-as.vector(EC[[1]]))))
  if (is.vector(z1)) {as.matrix(z1); Z <-t(sm12%*%z1) }
  else {Z<-t(sm12%*%t(z1))}

  return(list(y=Z, cumY = ECz))
}



#############################################*

# if true cum are known use to standardize
.Coef_GC<-function(cum){
  d <- length(cum[[1]])
  EC<-cum
  ECz <- vector(mode="list", length=4)
  ECz[[1]]=rep(0,d)
  ECz[[2]]=diag(d)
  #z1<-t(apply(X,1, function(x) x-as.vector(EC[[1]])))
  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx <- svd(cx)
  sm12 <- svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  #if (is.vector(z1)) {as.matrix(z1);Z<-t(sm12%*%z1) }  else {Z<-t(sm12%*%t(z1))}
  kron <- kronecker(sm12,sm12)

  for (j in c(3:4)){
    kron <- kronecker(kron,sm12)
    ECz[[j]]<- kron%*%EC[[j]]
  }
  return(ECz)
}

#############################################

.Coef_Ed<-function(cum,n){
  d <- length(cum[[1]])
  EC<-cum
  ECz <- vector(mode="list", length=4)
  ECz[[1]]=rep(0,d)
  ECz[[2]]=diag(d)
  #z1<-t(apply(X,1, function(x) x-as.vector(EC[[1]])))
  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx <- svd(cx)
  sm12 <- svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  #if (is.vector(z1)) {as.matrix(z1);Z<-t(sm12%*%z1) }  else {Z<-t(sm12%*%t(z1))}
  kron <- kronecker(sm12,sm12)

  for (j in c(3:4)){
    kron <- kronecker(kron,sm12)
    ECz[[j]]<- kron%*%EC[[j]]
  }

  CoeffEd <- vector(mode="list", length=3)
  CoeffEd[[1]] <- n^(-1/2) * ECz[[3]] / 6
  CoeffEd[[2]] <- n^(-1) * ECz[[4]] / 24
  CoeffEd[[3]] <- n^(-1) * kronecker(ECz[[3]],ECz[[3]]) / 72

  return(CoeffEd)
}

#############################################*






#########################################################
## Rcpp Fast functions

#' IntHermiteN
#'
#' Computes the integrals of d-Hermite polynomial with respect to the normal density
#' \deqn{\int_{-\infty }^{\mathbf{y}}\mathbf{H}_{k-1}\left( \mathbf{s}\right) \varphi
#' \left( \mathbf{s}\right) d \,\mathbf{s}}
#' either from \eqn{-\infty} to \eqn{\mathbf{y}} (type = "lower")
#' or from \eqn{\mathbf{y}} to \eqn{+\infty} ( type = "upper").
#'
#' @param x Numeric vector of lenght \eqn{d}
#' @param K Integer. The order of the Hermite polynomial + 1.
#' @param type Character string specifying the integration range.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"lower"}: integrate from \eqn{-\infty} to \eqn{x}
#'     \item \code{"upper"}: integrate from \eqn{x} to \eqn{+\infty}
#'   }
#'
#'
#' @return A list of integrated Hermite polynomials up to order K-1.
#'
#' @examples
#' x <- c(1,2)
#' IntHermiteN(x, K = 3,type = "lower")
#'
#' @export
IntHermiteN <- function(x, K, type = c("lower", "upper")) {
  type <- match.arg(type)

  # Basic input checks
  stopifnot(is.numeric(x), length(K) == 1, is.numeric(K))

  if (type == "lower") {
    res <- .Funct_CornF_CP(x, K)        # integrate from -Inf to x
  } else {
    res <- .Funct_CornF_low_CP(x, K)    #  integrate from x to +Inf
  }

  return(res)
}






## lower part see hermite_cf.ccp in \src
.Funct_CornF_CP <- function(x, K) {


  if (is.numeric(x) && is.null(dim(x))) {
    d <- length(x)
  } else {
    stop("x should be a numeric vector of length d")
  }

  mu   <- 0
  sig2 <- 1

  Fi.x <- pnorm(x, mean = mu, sd = sqrt(sig2))
  fi.x <- dnorm(x, mean = mu, sd = sqrt(sig2))

  # Hfi[k, m] = hermite(x[m], k) * fi.x[m], k = 1..K-2
  if (K > 2L) {
    Hfi <- matrix(0, nrow = K - 2L, ncol = d)
    for (k in 1:(K - 2L)) {
      for (m in 1:d) {
        Hfi[k, m] <- EQL::hermite(x[m], k) * fi.x[m]
      }
    }
  } else {
    Hfi <- matrix(0, nrow = 0L, ncol = d)
  }

  IntHermiteN_stream_cpp(Fi.x, fi.x, Hfi, K)
}

## upper part see IntHermiteN_stream_upper.ccp in \src
.Funct_CornF_low_CP <- function(x, K) {


  if (is.numeric(x) && is.null(dim(x))) {
    d <- length(x)
  } else {
    stop("x must be a numeric vector")
  }

  mu   <- 0
  sig2 <- 1

  Fi.x <- pnorm(x, mean = mu, sd = sqrt(sig2))
  fi.x <- dnorm(x, mean = mu, sd = sqrt(sig2))

  # Hfi[k, m]
  if (K > 2L) {
    Hfi <- matrix(0, nrow = K - 2L, ncol = d)
    for (k in 1:(K - 2L)) {
      for (m in 1:d) {
        Hfi[k, m] <- EQL::hermite(x[m], k) * fi.x[m]
      }
    }
  } else {
    Hfi <- matrix(0, nrow = 0L, ncol = d)
  }

  IntHermiteN_stream_upper_cpp(Fi.x, fi.x, Hfi, K)
}


#' Integrate Gram Charlier density
#'
#' Computes the integrals of the \eqn{d}-variate Gram Charlier  density
#' with respect to the normal density
#' It integrates the \code{\link{GramCharlier}} with the first 4
#' cumulants.
#'
#' @param x An nxd data matrix
#' @param cum Unstandardized first four cumlants
#' @param type Character string specifying the integration range.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"lower"}: integrate from \eqn{-\infty} to \eqn{x}
#'     \item \code{"upper"}: integrate from \eqn{x} to \eqn{+\infty}
#'   }
#'
#'
#' @return The vector of evaluated probabilities
#'
#'
#' @examples
#' x <- matrix(1:6,2,3, byrow=TRUE)
#' cum <- MomCumMVt(p = 12, d = 3, r = 4, nCum = TRUE)
#' # P(X <= x)
#' p <- IntGramCharlier(x, cum, type = "lower")
#'
#' @family Approximations
#'
#' @export
IntGramCharlier <- function(x, cum , type = c("lower", "upper")) {
  type <- match.arg(type)


  if (type == "lower") {
    res <- .CornishFisherM_CP(x,cum)        # integrate from -Inf to x
  } else {
    res <- .CornishFisherM_low_CP(x,cum)    # integrate from x to +Inf
  }

  return(res)
}

##################################################

#' Integrate Edgeworth density
#'
#' Computes the integrals of the \eqn{d}-variate Edgeworth  density
#' with respect to the normal density
#' It integrates the \code{\link{Edgeworth}} with the first 4
#' cumulants.
#'
#' @param x An nxd data matrix
#' @param cum Unstandardized first four cumlants
#' @param n the number of elements in the sum of data
#' @param type Character string specifying the integration range.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"lower"}: integrate from \eqn{-\infty} to \eqn{x}
#'     \item \code{"upper"}: integrate from \eqn{x} to \eqn{+\infty}
#'   }
#'
#'
#' @return  The vector of evaluated probabilities
#'
#'
#' @examples
#' x <- matrix(1:6,2,3, byrow=TRUE)
#' cum <- MomCumMVt(p = 12, d = 3, r = 4, nCum = TRUE)
#' # P(X <= x)
#' p <- IntEdgeworth(x, cum, type = "lower")
#'
#' @family Approximations
#'
#' @export
IntEdgeworth <- function(x, cum , n = 1, type = c("lower", "upper")) {
  type <- match.arg(type)


  if (type == "lower") {
    res <- .CorFisEdgeM_CP (x, cum, n)        # integrate from -Inf to x
  } else {
    res <- .CorFisEdgeM_low_CP (x, cum, n)    # integrate from x to +Inf
  }

  return(res)
}
########################################



.CornishFisherM_CP <- function(x,cum) {
  #   #,type
  #   #* input: x d-dimensional points for estimated distribution
  #   #* input: cum cumulants  the first 4 - 6 cumulants depends on type
  #   #* input: type "GC" or " Edgworth"
  #   #*  x <- Xt
  #   # cum <- cumt
  #Y.cumY  <- NULL
  Y.cumY  <- .X2Ystand(x,cum)
  Y <- Y.cumY[[1]]
  cumY <- Y.cumY[[2]]

  d <-  dim(x)[2]

  PHI <- apply(x, 1, function(x_row) {       # y_row <- x[1,]
    mvtnorm::pmvnorm(upper = x_row, mean = cum[[1]],
                     sigma = matrix(cum[[2]],nrow=d))[1]})

  pX <- apply(Y, 1, function(y_row){.CornishFisherM_in_CP(y_row,cumY)})+ PHI
  return(pX)
}

#############################################*

.CornishFisherM_low_CP <- function(x,cum) {
  #   #,type
  #   #* input: x d-dimensional points for estimated distribution
  #   #* input: cum cumulants  the first 4 - 6 cumulants depends on type
  #   #* input: type "GC" or " Edgworth"
  #   #*  x <- Xt
  #   # cum <- cumt
  #Y.cumY  <- NULL
  Y.cumY  <- .X2Ystand(x,cum)
  Y <- Y.cumY[[1]]
  cumY <- Y.cumY[[2]]

  d <-  dim(x)[2]

  PHI <- apply(x, 1, function(x_row) {       # y_row <- x[1,]
    mvtnorm::pmvnorm(lower = x_row, mean = cum[[1]],
                     sigma = matrix(cum[[2]],nrow=d))[1]})

  pX <- apply(Y, 1, function(y_row){.CornishFisherM_in_low_CP(y_row,cumY)})+ PHI
  return(pX)
}


#############################################*

.CorFisEdgeM_CP <- function(x,cum,n) {

  #Y.cumY  <- NULL
  Y.cumY  <- .X2Ystand(x,cum)
  Y <- Y.cumY[[1]]
  cumY <- Y.cumY[[2]]

  d <-  dim(x)[2]

  PHI <- apply(x, 1, function(x_row) {       # y_row <- x[1,]
    mvtnorm::pmvnorm(upper = x_row, mean = cum[[1]],
                     sigma = matrix(cum[[2]],nrow=d))[1]})

  pX <- apply(Y, 1, function(y_row){.CorFisEdgeM_in_CP(y_row,cumY,n)})+ PHI
  return(pX)
}

#############################################*


.CorFisEdgeM_low_CP <- function(x, cum, n) {
  #
  Y.cumY  <- .X2Ystand(x,cum)
  Y <- Y.cumY[[1]]
  cumY <- Y.cumY[[2]]

  d <-  dim(x)[2]

  PHI <- apply(x, 1, function(x_row) {       # y_row <- x[1,]
    mvtnorm::pmvnorm(lower = x_row, mean = cum[[1]],
                     sigma = matrix(cum[[2]],nrow=d))[1]})

  pX <- apply(Y, 1, function(y_row){.CorFisEdgeM_in_low_CP(y_row, cumY, n)})+ PHI
  return(pX)
}




#############################################*



.CornishFisherM_in_CP <- function(y_row,cumY) {
  K <- length(cumY)+1
  CF.SiH <- .Funct_CornF_CP(y_row,K)
  coeff.GC <- .Coef_GC(cumY)
  dGC <-   CF.SiH[[4]]%*%coeff.GC[[3]]/6 +
    CF.SiH[[5]]%*%coeff.GC[[4]]/24
  return(dGC)

}


#############################################*
.CornishFisherM_in_low_CP <- function(y_row,cumY) {
  K <- length(cumY)+1
  CF.SiH <- .Funct_CornF_low_CP(y_row,K)
  coeff.GC <- .Coef_GC(cumY)
  dGC <-   CF.SiH[[4]]%*%coeff.GC[[3]]/6 +
    CF.SiH[[5]]%*%coeff.GC[[4]]/24
  return(dGC)

}

#############################################*

.CorFisEdgeM_in_CP <- function(y_row,cumY,n) {
  K <- 7
  CF.SiH <- .Funct_CornF_CP(y_row,K)
  coeff.Ed <- .Coef_Ed(cumY,n)
  dEd <-   CF.SiH[[4]]%*%coeff.Ed[[1]] +
    CF.SiH[[5]]%*%coeff.Ed[[2]] +
    CF.SiH[[7]]%*%coeff.Ed[[3]]
  return(dEd)

}

#############################################*

.CorFisEdgeM_in_low_CP <- function(y_row,cumY,n) {
  K <- 7
  CF.SiH <- .Funct_CornF_low_CP(y_row,K)
  coeff.Ed <- .Coef_Ed(cumY,n)
  dEd <-   CF.SiH[[4]]%*%coeff.Ed[[1]] +
    CF.SiH[[5]]%*%coeff.Ed[[2]] +
    CF.SiH[[7]]%*%coeff.Ed[[3]]
  return(dEd)

}


#' Multivariate tail conditional expectation
#'
#' It provides the conditional expectation
#' \deqn{ \text{MTCE}_q(\mathbf{X}) = \operatorname{E} \left( \mathbf{X} \mid X_1 > \text{VaR}_q (X_1),
#'                           X_2 > \text{VaR}_q (X_2), \dots, X_n > \text{VaR}_q (X_d) \right),}
#' for \eqn{q \in (0,1)}, where \eqn{\text{VaR}_q(X)} is the q-th quantile of the random variable \eqn{X}.
#' Expectation is taken with respect to \code{\link{GramCharlier}} with the first 4
#' cumulants.
#'
#' For further details see the references below,
#'
#'
#' @param X a vector of unstandardized VaRq
#' @param cum list of mean, variance, skewness and kurtosis vectors
#'
#' @return Numerator of the ratio
#' @return Denominator of the ratio
#' @return MTCE Conditional expected value
#'
#' @references Landsman, Z., Makov, U., & Shushi, T. (2016).
#' Multivariate tail conditional expectation for elliptical
#' distributions. Insurance: Mathematics and
#' Economics, 70, 216-223.
#'
#' @examples
#' x <- c(2,3,4)
#' cum <- MomCumMVt(p = 12, d = 3, r = 4, nCum = TRUE)
#' CE <- MTCE(x, cum)
#'
#' @family Approximations
#'
#' @export
MTCE <- function (X, cum)
{
  k <- min(4, length(cum))
  if (k < 3)
    stop("k must be greater than 2")
  if (is.matrix(X))
    stop(" X must be a vector of VaRq values")
  d <- length(X)
  EC <- cum
  mu <- EC[[1]]
  ECz <- vector(mode = "list", length = k)
  ECz[[1]] = rep(0, d)
  ECz[[2]] = diag(d)
  z1 <- X - as.vector(EC[[1]])
  if (is.vector(EC[[2]])) {
    cx <- matrix(EC[[2]], nrow = d)
  }
  else {
    cx <- EC[[2]]
  }
  svdx <- svd(cx)
  sm12 <- svdx$u %*% diag(1/sqrt(svdx$d)) %*% t(svdx$u)
  s12 <- svdx$u %*% diag(sqrt(svdx$d)) %*% t(svdx$u)
  Z <- c(sm12 %*% z1)

  IntHm_U <- .Funct_CornF_low_CP(Z,6)

  Int1 <- NULL
  for (i in 1:d) {
    Int1[i] <- dnorm(Z[i]) * mvtnorm::pmvnorm(lower = Z[-i])
  }
  ST1 <- IntHm_U[[5]] + 3 * (kronecker(SymMatr(d, 3), diag(d))) %*%
    kronecker(IntHm_U[[3]], c(diag(d)))
  M3 <- matrix(ST1, d, d^3)
  SK <- c(M3 %*% EC[[3]])
  KT1 <- IntHm_U[[6]] + 4 * (kronecker(SymMatr(d, 4), diag(d))) %*%
    kronecker(IntHm_U[[4]], c(diag(d)))
  M4 <- matrix(KT1, d, d^4)
  KK <- c(M4 %*% EC[[4]])
  Num <- Int1 + SK + KK
  Int2 <- mvtnorm::pmvnorm(lower = Z)
  DSK <- c(IntHm_U[[4]] %*% EC[[3]])
  DKK <- c(IntHm_U[[5]] %*% EC[[4]])
  Den <- Int2 + DSK + DKK
  lbdq <- Num/Den
  nMTCE <- c(mu + s12 %*% lbdq)
  output <- list(Numerator = Num, Denominator = Den, Lambda = lbdq,
                 MTCE = nMTCE)
  return(output)
}


