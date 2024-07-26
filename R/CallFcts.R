# 1. CommutatorMatr
# 2. Partitions
# 3. CommutatorIndx
# 4. HermiteN2X
# 5. HermiteN
# 6. HermiteCoeff
# 7. Mom2Cum
# 8. Cum2Mom
# 9. SampleKurtosis
#10. SampleSkew
#11. MomCumZabs
#12. EVSKUniS



#' Commutator Matrix
#'
#' This function generates various types of commutator matrices.
#'
#' @param Type A string specifying the type of commutator matrix. Choices are "Kmn", "Kperm", "Mixing", or "Moment".
#' @param ... Additional arguments specific to the type of commutator matrix (see Details).
#'
#' @details
#' The function \code{CommutatorMatr} supports the following types of commutator matrices:
#'
#' \describe{
#'   \item{Kmn}{
#'     \strong{Description:} Transforms \code{vec(A)} to \code{vec(A^T)}, where \code{A^T} is the transpose of matrix \code{A}. An option for sparse matrix is provided. By default, a non-sparse matrix is produced. Using sparse matrices increases computation times but requires far less memory.
#'     \strong{Arguments:}
#'     \describe{
#'       \item{\code{m} (integer)}{Number of rows of the first matrix.}
#'       \item{\code{n} (integer)}{Number of columns of the first matrix.}
#'       \item{\code{useSparse} (logical, optional)}{If TRUE, returns a sparse matrix. Default is FALSE.}
#'     }
#'   }
#'   \item{Kperm}{
#'     \strong{Description:} Generates a commutation matrix for a specified permutation of matrix dimensions. An option for sparse matrix is provided. By default, a non-sparse matrix is produced. Using sparse matrices increases computation times but requires far less memory.
#'     \strong{Arguments:}
#'     \describe{
#'       \item{\code{perm} (integer vector)}{The permutation vector.}
#'       \item{\code{dims} (integer vector)}{The dimensions of the matrices involved.}
#'       \item{\code{useSparse} (logical, optional)}{If TRUE, returns a sparse matrix. Default is FALSE.}
#'     }
#'   }
#'   \item{Mixing}{
#'     \strong{Description:} Generates the Mixing commutation matrix used in linear algebra transformations involving tensor products. An option for sparse matrix is provided. By default, a non-sparse matrix is produced. Using sparse matrices increases computation times but requires far less memory.
#'     \strong{Arguments:}
#'     \describe{
#'       \item{\code{d1} (integer vector)}{Dimensions of the first set.}
#'       \item{\code{d2} (integer vector)}{Dimensions of the second set.}
#'       \item{\code{useSparse} (logical, optional)}{If TRUE, returns a sparse matrix. Default is FALSE.}
#'     }
#'   }
#'   \item{Moment}{
#'     \strong{Description:} Generates the Moment commutation matrix based on partitioning of moments. An option for sparse matrix is provided. By default, a non-sparse matrix is produced. Using sparse matrices increases computation times but requires far less memory.
#'     \strong{Arguments:}
#'     \describe{
#'       \item{\code{el_rm} (integer vector)}{Elements of the partition.}
#'       \item{\code{d} (integer)}{Dimension of the partition.}
#'       \item{\code{useSparse} (logical, optional)}{If TRUE, returns a sparse matrix. Default is FALSE.}
#'     }
#'   }
#' }
#'
#' @return Depending on the  type:
#' \describe{
#'   \item{Kmn}{A commutation matrix of dimension \eqn{mn \times mn}. If \code{useSparse=TRUE}, an object of class "dgCMatrix" is produced.}
#'   \item{Kperm}{A square permutation matrix of size \code{prod(dims)}. If \code{useSparse=TRUE}, an object of class "dgCMatrix" is produced.}
#'   \item{Mixing}{A square matrix of dimension \code{prod(d1) * prod(d2)}. If \code{useSparse=TRUE}, an object of class "dgCMatrix" is produced.}
#'   \item{Moment}{A commutator matrix for moment formulae.}
#' }
#' @examples
#' # Example for Kmn
#' CommutatorMatr("Kmn", m = 3, n = 2)
#'
#' # Example for Kperm
#' dims <- c(2, 3, 2)
#' perm <- c(1, 3, 2)
#' CommutatorMatr("Kperm", perm = perm, dims = dims)
#'
#' # Example for Mixing
#' d1 <- c(2, 3, 2)
#' d2 <- c(3, 2, 2)
#' CommutatorMatr("Mixing", d1 = d1, d2 = d2)
#'
#' # Example for Moment
#' n <- 4
#' r <- 2
#' m <- 1
#' d <- 2
#' PTA <- PartitionTypeAll(n)
#' el_r <- PTA$eL_r[[r]][m,]
#' CommutatorMatr("Moment", el_r = el_r, d = d)
#' @family Commutators
#' @export

CommutatorMatr <- function(Type, ...) {
  if (Type == "Kmn") {
    return(.matr_Commutator_Kmn(...))
  } else if (Type == "Kperm") {
    return(.matr_Commutator_Kperm(...))
  } else if (Type == "Mixing") {
    return(.matr_Commutator_Mixing(...))
  } else if (Type == "Moment") {
    return(.matr_Commutator_Moment(...))
  } else {
    stop("Invalid type. Choose from 'Kmn', 'Kperm', 'Mixing', 'Moment'.")
  }
}



#' General Partition Function
#'
#' A unified function to compute different types of partitions. Depending on the partition type specified, it calls the appropriate function: Partition_2Perm, Partition_DiagramsClosedNoLoops, Partition_Indecomposable, or Partition_Pairs.
#'
#' @param Type A character string specifying the type of partion to compute. Choose from "2Perm", "Diagram", "Indecomp", "Pairs".
#' @param ... Additional arguments passed to the specific partition function:
#' \describe{
#'   \item{{For "2Perm", "Diagram" and "Indecomp":}}{
#'     \itemize{
#'       \item{\code{L}}: A partition matrix.
#'     }
#'   }
#'   \item{{For "Pairs":}}{
#'     \itemize{
#'       \item{\code{N}}: An integer specifying the number of elements to be partitioned.
#'     }
#'   }
#' }
#' @return Depending on the commutator type:
#' \describe{
#'   \item{2Perm}{A vector with the elements 1 to N permuted according to L.}
#'   \item{Diagram}{The list of partition matrices indecomposable with respect to L, representing diagrams without loops.}
#'   \item{Indecomp}{A list of partition matrices indecomposable with respect to L and a vector indicating the number of indecomposable partitions by sizes.}
#'   \item{Pairs}{The list of partition matrices with blocks containing two elements. The list is empty if N is odd.}
#' }
#' @examples
#' # Example for 2Perm
#' PA <- PartitionTypeAll(4)
#' Partitions("2Perm", L = PA$Part.class[[3]])
#'
#' # Example for Diagram
#' L <- matrix(c(1,1,0,0,0,0,1,1),2,4,byrow=TRUE)
#' Partitions("Diagram", L = L)
#'
#' # Example for Indecomp
#' L <- matrix(c(1,1,0,0,0,0,1,1),2,4,byrow=TRUE)
#' Partitions("Indecomp", L = L)
#'
#' # Example for Pairs
#' Partitions("Pairs", N = 4)
#' @family Partitions
#' @export
Partitions <- function(Type, ...) {
  args <- list(...)
  if (Type == "2Perm") {
    if (!"L" %in% names(args)) stop("Argument 'L' is required for partition type '2Perm'")
    return(.Partition_2Perm(args$L))
  } else if (Type == "Diagram") {
    if (!"L" %in% names(args)) stop("Argument 'L' is required for partition type 'Diagram'")
    return(.Partition_DiagramsClosedNoLoops(args$L))
  } else if (Type == "Indecomp") {
    if (!"L" %in% names(args)) stop("Argument 'L' is required for partition type 'Indecomp'")
    return(.Partition_Indecomposable(args$L))
  } else if (Type == "Pairs") {
    if (!"N" %in% names(args)) stop("Argument 'N' is required for partition type 'Pairs'")
    return(.Partition_Pairs(args$N))
  } else {
    stop("Invalid partition type specified.")
  }
}



#' Commutator Index
#'
#' This function calculates the commutator index based on the specified type. The available types are "Kmn", "Kperm", "Mixing", and "Moment".
#' Depending on the selected type, the corresponding specific function is called.
#'
#' @param Type a string specifying the type of commutator index to be calculated. Must be one of "Kmn", "Kperm", "Mixing", or "Moment".
#' @param ... additional arguments passed to the specific commutator function.
#'
#' @return A vector representing the commutator index.
#'
#' @details
#' The function `CommutatorIndx` acts as a wrapper to call specific commutator functions based on the input `Type`.
#'
#' \strong{Type "Kmn"}:
#' \describe{
#'   \item{Parameters:}{
#'     \itemize{
#'       \item \code{m} - Row-dimension.
#'       \item \code{n} - Col-dimension.
#'     }
#'   }
#'   \item{Return:}{
#'     A vector of indexes to provide the commutation, transforming vec A to vec of the transposed A.
#'   }
#' }
#'
#' \strong{Type "Kperm"}:
#' \describe{
#'   \item{Parameters:}{
#'     \itemize{
#'       \item \code{perm} - Vector indicating the permutation of the order in the Kronecker product.
#'       \item \code{dims} - Vector indicating the dimensions of the vectors.
#'     }
#'   }
#'   \item{Return:}{
#'     An index vector to produce the permutation of the Kronecker products of vectors of any length.
#'   }
#' }
#'
#' \strong{Type "Mixing"}:
#' \describe{
#'   \item{Parameters:}{
#'     \itemize{
#'       \item \code{x} - A vector of dimension \code{prod(d1)*prod(d2)}.
#'       \item \code{d1} - Dimension of the first group of vectors.
#'       \item \code{d2} - Dimension of the second group of vectors.
#'     }
#'   }
#'   \item{Return:}{
#'     A vector Kx representing the product of the moment commutator and the vector x.
#'   }
#' }
#'
#' \strong{Type "Moment"}:
#' \describe{
#'   \item{Parameters:}{
#'     \itemize{
#'       \item \code{x} - A vector of length \code{d^n} where n is the length of \code{el_rm}.
#'       \item \code{el_rm} - Type of a partition.
#'       \item \code{d} - Dimensionality of the underlying multivariate distribution.
#'     }
#'   }
#'   \item{Return:}{
#'     A vector Kx representing the product of the moment commutator and the vector x.
#'   }
#' }
#'
#' @examples
#' # Kmn example
#' A <- 1:6
#' CommutatorIndx(Type = "Kmn", m = 3, n = 2)
#'
#' # Kperm example
#' a1 <- c(1, 2)
#' a2 <- c(2, 3, 4)
#' a3 <- c(1, 3)
#' p1 <- a1 %x% a2 %x% a3
#' CommutatorIndx(Type = "Kperm", perm = c(3, 1, 2), dims = c(2, 3, 2))
#'
#' # Mixing example
#' d1 <- c(2, 3, 2)
#' d2 <- c(3, 2, 2)
#' x <- 1:(prod(d1) * prod(d2))
#' CommutatorIndx(Type = "Mixing", x = x, d1 = d1, d2 = d2)
#'
#' # Moment example
#' n <- 4
#' r <- 2
#' m <- 1
#' d <- 2
#' PTA <- PartitionTypeAll(n)
#' el_r <- PTA$eL_r[[r]][m, ]
#' x <- 1:d^n
#' CommutatorIndx(Type = "Moment", x = x, el_rm = el_r, d = d)
#'
#' @family Commutators
#' @export
CommutatorIndx <- function(Type, ...) {
  switch(Type,
         "Kmn" = .indx_Commutator_Kmn(...),
         "Kperm" = .indx_Commutator_Kperm(...),
         "Mixing" = .indx_Commutator_Mixing(...),
         "Moment" = .indx_Commutator_Moment(...),
         stop("Invalid Type. Choose from 'Kmn', 'Kperm', 'Mixing', 'Moment'."))
}





########################################

#' Inverse Hermite Polynomial
#'
#' Compute the inverse of univariate or multivariate Hermite polynomials.
#'
#' This function computes the powers of x when Hermite polynomials are given.
#' Depending on the type specified, it handles either univariate or multivariate
#' Hermite polynomials.
#'
#' @param Type A string specifying the type of Hermite polynomial inversion.
#' Must be either "Univariate" or "Multivariate".
#' @param H_N Input Hermite polynomials. For univariate, it is a vector. For multivariate, it is a list.
#' @param N The highest polynomial order.
#' @param Sig2 The variance matrix of x for multivariate, or variance for univariate. Defaults to identity matrix for multivariate and 1 for univariate.
#'
#' @return A list of x powers: \eqn{x}, \eqn{x^{\otimes 2}}, ... , \eqn{x^{\otimes N}} for multivariate,
#' or a vector of x powers: \eqn{x^n}, \eqn{n=1:N} for univariate.
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. Section 4.6.2, (4.72), p.223 and Section 4.4, (4.23), p.198.
#'
#' @examples
#' # Univariate example
#' H_N_x <- c(1, 2, 3, 4)
#' x_powers <- HermiteN2X(Type = "Univariate", H_N = H_N_x, N = 4, Sig2 = 1)
#'
#' # Multivariate example
#' x <- c(1, 3)
#' Sig2 <- diag(length(x))
#' N <- 4
#' H_N_X <- HermiteN(x, N, Type="Multivariate")
#' x_ad_n <- HermiteN2X(Type = "Multivariate", H_N = H_N_X, N = N, Sig2 = Sig2)
#'
#' @family Hermite Polynomials
#' @export

HermiteN2X <- function(Type, H_N, N, Sig2 = NULL) {
  if (Type == "Univariate") {
    if (is.null(Sig2)) {
      Sig2 <- 1
    }
    return(.Hermite_Poly_NH_Inv(H_N, Sig2))
  } else if (Type == "Multivariate") {
    if (is.null(Sig2)) {
      Sig2 <- diag(length(H_N[[1]]))
    }
    return(.Hermite_Poly_NH_Multi_Inv(H_N, N, Sig2))
  } else {
    stop("Invalid Type. Must be either 'Univariate' or 'Multivariate'.")
  }
}

#

#' Hermite Polynomials (Univariate and Multivariate)
#'
#' Computes either univariate or multivariate Hermite polynomials up to a specified order.
#'
#' Depending on the value of the `Type` parameter, this function computes either the univariate or the multivariate Hermite polynomials.
#'
#' @param x A scalar (for univariate) or a vector (for multivariate) at which to evaluate the Hermite polynomials.
#' @param N The maximum order of the polynomials.
#' @param Type A character string specifying the type of Hermite polynomials to compute. Can be either "Univariate" or "Multivariate".
#' @param sigma2 The variance for univariate Hermite polynomials. Default is 1. (Only used if Type is "Univariate").
#' @param Sig2 The covariance matrix for multivariate Hermite polynomials. Default is the unit matrix diag(length(x)). (Only used if Type is "Multivariate").
#'
#' @return Depending on the type, the function returns:
#' \itemize{
#'   \item \code{Univariate}: A vector of univariate Hermite polynomials with degrees from 1 to N evaluated at x.
#'   \item \code{Multivariate}: A list of multivariate polynomials of order from 1 to N evaluated at vector x.
#' }
#'
#' @references
#' Gy.Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Section 4.1 (for univariate), Section 4.6.2, (4.73), p.223 (for multivariate).
#'
#' @examples
#' # Univariate example
#' HermiteN(x = 1, N = 3, Type = "Univariate")
#'
#' # Multivariate example
#' HermiteN(x = c(1, 3), N = 3, Type = "Multivariate", Sig2 = diag(2))
#'
#' @family Hermite Polynomials
#'
#' @export
HermiteN <- function(x, N, Type, sigma2 = 1, Sig2 = diag(length(x))) {
  if (Type == "Univariate") {
    return(.Hermite_Poly_HN(x, N, sigma2))
  } else if (Type == "Multivariate") {
    return(.Hermite_Poly_HN_Multi(x, N, Sig2))
  } else {
    stop("Invalid Type. Please specify either 'Univariate' or 'Multivariate'.")
  }
}



#' Coefficients of Hermite polynomials
#'
#' Provides the coefficients of Hermite polynomials, either univariate or multivariate.
#'
#' @param Type A character string specifying the type of Hermite polynomial. Must be either "Univariate" or "Multivariate".
#' @param N The order of polynomial. Required for both types.
#' @param d The dimension of the d-variate X. Required only for multivariate type.
#' @return For `Type = "Univariate"`, returns a vector of coefficients of \eqn{x^N}, \eqn{x^{N-2}}, etc.
#' For `Type = "Multivariate"`, returns a list of matrices of coefficients for the d-variate polynomials from 1 to N.
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Sections 4.4 (4.24) and 4.6.2, p. 223, Remark 4.8
#'
#' @family Hermite Polynomials
#' @export
#'
#' @examples
#' # Univariate example
#' H_uni <- HermiteCoeff(Type = "Univariate", N = 5)
#'
#' # Multivariate example
#' N <- 5; d <- 3
#' H_multi <- HermiteCoeff(Type = "Multivariate", N = N, d = d)
#' X <- c(1:3)
#' X3 <- kronecker(X, kronecker(X, X))
#' X5 <- kronecker(X3, kronecker(X, X))
#' Idv <- as.vector(diag(d)) # vector of variance matrix
#' # value of H5 at X is
#' vH5 <- H_multi[[1]] %*% X5 + H_multi[[2]] %*% kronecker(Idv, X3) +
#'   H_multi[[3]] %*% kronecker(kronecker(Idv, Idv), X)
HermiteCoeff <- function(Type, N, d = NULL) {
  if (Type == "Univariate") {
    return(.Hermite_Coeff(N))
  } else if (Type == "Multivariate") {
    if (is.null(d)) {
      stop("Parameter 'd' must be provided for multivariate Hermite coefficients")
    }
    return(.Hermite_CoeffMulti(N, d))
  } else {
    stop("Invalid Type. Must be either 'Univariate' or 'Multivariate'.")
  }
}


#' Convert moments to cumulants (univariate and multivariate)
#'
#' Obtains a vector of cumulants from a vector of moments for either univariate or multivariate data.
#'
#' @param moments Either a vector of univariate moments or a list of vectors of multivariate moments.
#' @param Type A character string specifying the type of moments provided. Use "Univariate" for univariate moments and "Multivariate" for multivariate moments.
#'
#' @return The vector of cumulants if \code{Type} is "Univariate" or the list of vectors of cumulants if \code{Type} is "Multivariate".
#'
#' @examples
#' # Univariate example
#' mu_x <- c(1, 2, 3, 4)
#' Mom2Cum(mu_x, Type = "Univariate")
#'
#' # Multivariate example
#' mu <- list(c(0,0), c(1,0,0,1), c(0,0,0,0,0,0,0,0), c(3,0,0,1,0,1,1,0,0,1,1,0,1,0,0,3), c(rep(0,32)))
#' Mom2Cum(mu, Type = "Multivariate")
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Section 3.4.
#'
#' @family Moments and cumulants
#'
#' @export
Mom2Cum <- function(moments, Type = c("Univariate", "Multivariate")) {
  Type <- match.arg(Type)

  if (Type == "Univariate") {
    .conv_Mom2Cum(moments)
  } else if (Type == "Multivariate") {
    .conv_Mom2CumMulti(moments)
  } else {
    stop("Invalid Type. Use 'Univariate' or 'Multivariate'.")
  }
}





#' Convert cumulants to moments (univariate and multivariate)
#'
#' Obtains a vector of moments from a vector of cumulants for either univariate or multivariate data.
#'
#' @param cumulants Either a vector of univariate cumulants or a list of vectors of multivariate cumulants.
#' @param Type A character string specifying the type of cumulants provided. Use "Univariate" for univariate cumulants and "Multivariate" for multivariate cumulants.
#'
#' @return The vector of moments if \code{Type} is "Univariate" or the list of vectors of moments if \code{Type} is "Multivariate".
#'
#' @examples
#' # Univariate example
#' cum_x <- c(1, 2, 3, 4)
#' Cum2Mom(cum_x, Type = "Univariate")
#'
#' # Multivariate example
#' cum <- list(c(0,0), c(1,0,0,1), c(rep(0,8)), c(rep(0,16)), c(rep(0,32)))
#' Cum2Mom(cum, Type = "Multivariate")
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Section 3.4.
#'
#' @family Moments and cumulants
#'
#' @export
Cum2Mom <- function(cumulants, Type = c("Univariate", "Multivariate")) {
  Type <- match.arg(Type)

  if (Type == "Univariate") {
    .conv_Cum2Mom(cumulants)
  } else if (Type == "Multivariate") {
    .conv_Cum2MomMulti(cumulants)
  } else {
    stop("Invalid Type. Use 'Univariate' or 'Multivariate'.")
  }
}


#' Estimation of Sample Kurtosis (Mardia, MRSz, Total)
#'
#' Estimates the sample kurtosis index based on the specified method: Mardia, MRSz, or Total.
#'
#' @param x A matrix of multivariate data.
#' @param Type A character string specifying the type of kurtosis index to estimate. Use "Mardia" for Mardia's kurtosis index, "MRSz" for the Mori-Rohatgi-Szekely kurtosis matrix, or "Total" for the total kurtosis index.
#'
#' @return A list containing the estimated kurtosis index or matrix and the associated p-value under the Gaussian hypothesis.
#' \item{Mardia.Kurtosis}{The kurtosis index when \code{Type} is "Mardia".}
#' \item{MRSz.Kurtosis}{The kurtosis matrix when \code{Type} is "MRSz".}
#' \item{Total.Kurtosis}{The total kurtosis index when \code{Type} is "Total".}
#' \item{p.value}{The p-value under the Gaussian hypothesis for the estimated kurtosis.}
#'
#' @examples
#' # Mardia's kurtosis example
#' x <- matrix(rnorm(100*5), ncol=5)
#' SampleKurt(x, Type = "Mardia")
#'
#' # MRSz's kurtosis example
#' SampleKurt(x, Type = "MRSz")
#'
#' # Total kurtosis example
#' SampleKurt(x, Type = "Total")
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Example 6.1 and 6.9.
#'
#' @family Estimation
#'
#' @export
SampleKurt <- function(x, Type = c("Mardia", "MRSz", "Total")) {
  Type <- match.arg(Type)

  if (Type == "Mardia") {
    .Esti_Kurt_Mardia(x)
  } else if (Type == "MRSz") {
    .Esti_Kurt_MRSz(x)
  } else if (Type == "Total") {
    .Esti_Kurt_Total(x)
  } else {
    stop("Invalid Type. Use 'Mardia', 'MRSz', or 'Total'.")
  }
}


#' Estimation of Sample Skewness (Mardia, MRSz)
#'
#' Estimates the sample skewness index based on the specified method: Mardia or MRSz.
#'
#' @param x A matrix of multivariate data.
#' @param Type A character string specifying the type of skewness index to estimate. Use "Mardia" for Mardia's skewness index or "MRSz" for the Mori-Rohatgi-Szekely skewness vector and index.
#'
#' @return A list containing the estimated skewness index or vector and the associated p-value under the Gaussian hypothesis.
#' \item{Mardia.Skewness}{The skewness index when \code{Type} is "Mardia".}
#' \item{MRSz.Skewness.Vector}{The skewness vector when \code{Type} is "MRSz".}
#' \item{MRSz.Skewness.Index}{The skewness index when \code{Type} is "MRSz".}
#' \item{p.value}{The p-value under the Gaussian hypothesis for the estimated skewness.}
#'
#' @examples
#' # Mardia's skewness example
#' x <- matrix(rnorm(100*5), ncol=5)
#' SampleSkew(x, Type = "Mardia")
#'
#' # MRSz's skewness example
#' SampleSkew(x, Type = "MRSz")
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear, Springer 2021. Example 6.1 and 6.2.
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate skewness and kurtosis. Sankhya A, 83(2), 607-644.
#'
#' @family Estimation
#'
#' @export
SampleSkew <- function(x, Type = c("Mardia", "MRSz")) {
  Type <- match.arg(Type)

  if (Type == "Mardia") {
    .Esti_Skew_Mardia(x)
  } else if (Type == "MRSz") {
    .Esti_Skew_MRSz(x)
  } else {
    stop("Invalid Type. Use 'Mardia' or 'MRSz'.")
  }
}


#' Moments and Cumulants of the Central Folded Normal Distribution
#'
#' Provides the theoretical moments and cumulants of the Central Folded Normal distribution.
#' Depending on the choice of `Type`, either the univariate or d-variate distribution is used.
#'
#' @param r The highest moment (cumulant) order.
#' @param nCum Logical; if TRUE, then cumulants are calculated.
#' @param Type Character; specifies the type of distribution. Must be either "Univariate" or "Multivariate".
#' @param d Integer; the dimension of the distribution. Must be 1 when `Type` is "Univariate" and greater than 1 when `Type` is "Multivariate".
#' @return A list containing moments and optionally cumulants.
#' \itemize{
#'   \item For "Univariate" type:
#'     \itemize{
#'       \item \code{MuZ}: The moments of the univariate Central Folded Normal distribution.
#'       \item \code{CumZ}: The cumulants of the univariate Central Folded Normal distribution.
#'     }
#'   \item For "Multivariate" type:
#'     \itemize{
#'       \item \code{MuZ}: The moments of the d-variate Central Folded Normal distribution.
#'       \item \code{CumZ}: The cumulants of the d-variate Central Folded Normal distribution.
#'     }
#' }
#' @examples
#' # Univariate case: The first three moments
#' MomCumZabs(3, 1, Type = "Univariate")
#' # Univariate case: The first three moments and cumulants
#' MomCumZabs(3, 1, Type = "Univariate",nCum = TRUE)
#' # d-variate case: The first three moments
#' MomCumZabs(3, 2, Type = "Multivariate" )
#' # d-variate case: The first three moments and cumulants
#' MomCumZabs(3, d=2, Type = "Multivariate",  nCum = TRUE)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Proposition 5.1 p.242 and formula: p. 301
#' @family Moments and cumulants
#' @export
#'
MomCumZabs <- function(r, d, Type, nCum = FALSE) {
  if (Type == "Univariate" && d != 1) {
    stop("For 'Univariate' Type, d must be 1.")
  } else if (Type == "Multivariate" && d <= 1) {
    stop("For 'Multivariate' Type, d must be greater than 1.")
  }

  if (Type == "Univariate") {
    return(.distr_Zabs_MomCum_Th(r, nCum))
  } else if (Type == "Multivariate") {
    return(.distr_ZabsM_MomCum_Th(r,d,nCum))
  } else {
    stop("Invalid Type. Please choose either 'Univariate' or 'Multivariate'.")
  }
}



#' EVSK of the Uniform distribution on the sphere or its modulus
#'
#' Cumulants (up to the 4th order), skewness, and kurtosis of the d-variate Uniform distribution on
#' the sphere or the modulus of the d-variate Uniform distribution on the sphere.
#'
#' @param d dimensions
#' @param nCum if it is FALSE then moments (up to the 4th order) are calculated.
#' @param Type specify the type of distribution: "Standard" for the Uniform distribution on the sphere,
#' or "Modulus" for the modulus of the Uniform distribution on the sphere.
#' @return A list of computed moments and cumulants.
#'
#' When Type is "Standard":
#' \item{EU1}{Mean vector}
#' \item{varU}{Covariance matrix}
#' \item{Skew.U}{Skewness vector (always zero)}
#' \item{Skew.tot}{Total skewness (always zero)}
#' \item{Kurt.U}{Kurtosis vector}
#' \item{Kurt.tot}{Total kurtosis}
#'
#' When Type is "Modulus":
#' \item{EU1}{Mean vector}
#' \item{varU}{Covariance matrix}
#' \item{EU.k}{List of moments up to 4th order}
#' \item{cumU.k}{List of cumulants up to 4th order}
#' \item{skew.U}{Skewness vector}
#' \item{kurt.U}{Kurtosis vector}
#'
#' @references Gy. Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 Proposition 5.3 p.297
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @references Gy. Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Lemma 5.12 p.298
#' @family Moments and cumulants
#' @export
#' @examples
#' # Example for Standard type
#' EVSKUniS(d=3, Type="Standard")
#'
#' # Example for Modulus type
#' EVSKUniS(d=3, Type="Modulus")
EVSKUniS <- function(d, nCum = TRUE, Type = c("Standard", "Modulus")) {
  Type <- match.arg(Type)

  if (Type == "Standard") {
    return(.distr_Uni_EVSK_Th(d, nCum))
  } else if (Type == "Modulus") {
    return(.distr_UniAbs_EVSK_Th(d, nCum))
  } else {
    stop("Invalid Type. Please choose either 'Standard' or 'Modulus'.")
  }
}
