###########
## 1. matr_Commutator_Kmn  1.2.3 Commutation Matrix p. 8, (1.12)
##    - With option for sparse matrices
## 2. matr_Commutator_Kperm 1.2.4 Commuting T-Products of Vectors p. 11, (1.23)
##  - With option for sparse matrices
## 3. matr_Commutator_Mixing; 4.6 Moments, Cumulants, and Linearization, p.218,(4.58)
##        A.2.2.1 Mixing Commutator
##  - With option for sparse matrices
## 4. Permutation_Inverse
## 5. SymMatr - calculates symmetrizer; 1.3.1 Symmetrization, p.14. (1.29)
## 6. EliminMatr 1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
## 7. QplicMatr, p.21, (1.31)
## 8. UnivMomCum - selects univariate moments and cumulants from the T-vector cumulants



### 1.2.3 Commutation Matrix p. 8, (1.12)





.matr_Commutator_Kmn <-  function(m,n,useSparse=FALSE){
  pp <- 1:(m*n)
  #p2 <- t(matrix(data=pp,nrow = m, ncol = n))
  p2 <- matrix(data=pp,nrow = n, ncol = m,byrow=TRUE)
  pc <- as.vector(p2)
  if (useSparse==FALSE) {Ic <- diag(m*n)}
  else
  {Ic <- Matrix::Diagonal(m*n)}
  Ic <- Ic[pc,]
  return(Ic)
}





.matr_Commutator_Kperm <- function(perm,dims,useSparse=FALSE){
  n <- length(perm)
  if (length(dims)== 1) {
    d<-dims
    #dims <- rep(dims,n)
    k <- length(perm)  # k <- length(permutation0)
    # egyk <- ones(k,1);
    perm0 <- perm    # perm0 <- permutation0
    perm01 <- rev(perm0)
    perm02 <-  order(perm01)
    permutation0 <- rev(perm02)

    pik <- function(perm,d) { 1+sum((perm-1)*cumprod(d*rep(1,k))/d)}  # pik számolás

    invIk0 <-  order(permutation0)#  invIk0  inverz permutació

    Allind <-  unique(arrangements::combinations( rep(1:d, k), k ) ) #= NULL, n = NULL
    # unique(nchoosek(repmat(1:d, 1,k), k), 'rows'); %all indeces
    if (useSparse==FALSE) {K_perm <- matrix(rep(0, (d^(2*k))), nrow= d^k)}
    else {K_perm <- Matrix::Matrix(rep(0, (d^(2*k))), nrow= d^k,sparse=TRUE)}
    for (jj in c(1:(d^k))){
      pA <- Allind[jj,]
      pInv <- pA[invIk0] #inverzepermutation0
      ind1 <- pik(pInv,d)
      ind2 <- pik(pA,d)
      K_perm[ind1,ind2] <- 1;
    }
    return(K_perm)}
  S.perm <-  sort(perm , index.return = TRUE)
  perm <- S.perm$ix # % inverse ( sort is increasing)

  if (useSparse==FALSE) {
    K.matr <- diag(prod(dims))
  if (length(perm)== 1) {return( K.matr )}
  U.before <- NULL
  U.after <- NULL
  for (i in 1 : (n - 1)) {
    # run loop (n-i) times
    for (j in 1 : (n - i)) {
      # compare elements
      if (perm[j] > perm[j + 1]) {
        if ((j-1)>0) {U.before <- diag(prod(dims[1:(j-1)]))}
        else {U.before <-1}

        if ((j+1)==n){ U.after <- 1 }
        else { U.after <- diag(prod(dims[(j+2):n]))}

        K1.matr <-  kronecker(U.before,kronecker(.matr_Commutator_Kmn(dims[j+1],dims[j]),U.after))
        # K2.matr <- (K1.matr,)
        K.matr <- K1.matr %*% K.matr
        temp.d <- dims[j]
        dims[j] <- dims[j+1]
        dims[j+1] <- temp.d
        temp.p <- perm[j]
        perm[j] <- perm[j + 1]
        perm[j + 1] <- temp.p
      }
    }
  }
  }

  if (useSparse==TRUE) {
    K.matr <- Matrix::Diagonal(prod(dims))
    if (length(perm)== 1) {return( K.matr )}
    U.before <- NULL
    U.after <- NULL
    for (i in 1 : (n - 1)) {
      # run loop (n-i) times
      for (j in 1 : (n - i)) {
        # compare elements
        if (perm[j] > perm[j + 1]) {
          if ((j-1)>0) {U.before <- Matrix::Diagonal(prod(dims[1:(j-1)]))}
          else {U.before <-1}

          if ((j+1)==n){ U.after <- 1 }
          else { U.after <- Matrix::Diagonal(prod(dims[(j+2):n]))}

          K1.matr <-  kronecker(U.before,kronecker(.matr_Commutator_Kmn(dims[j+1],dims[j],useSparse=TRUE),U.after))
          # K2.matr <- (K1.matr,)
          K.matr <- K1.matr %*% K.matr
          temp.d <- dims[j]
          dims[j] <- dims[j+1]
          dims[j+1] <- temp.d
          temp.p <- perm[j]
          perm[j] <- perm[j + 1]
          perm[j + 1] <- temp.p
        }
      }
    }
  }
  return(K.matr)
}







.matr_Commutator_Mixing <- function( d1,d2,useSparse=FALSE) {
  # permutations
  # dim(d1)<-dim(d2)
  n <- length(d2)
  if (useSparse==FALSE){
    i1<- matrix(data =c(1:n,(1:n)+n) , nrow = 2, ncol = n, byrow =TRUE)
    fact_n<-factorial(n)  # number of permutations
    B <- matrix(data = rep(c(1:n),fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
    Permut <- arrangements::permutations(1:n)#perm(c(1:n))
    q<-cbind(B, (Permut+n))
    indUj<- c(i1) # reorder q
    q<-q[,indUj] # permutations
    # dimensions
    d1B <- matrix(rep(d1,fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
    d2q<- matrix(rep(0,fact_n*n),nrow = fact_n, ncol = n) # zeros(fact_n,n);
      for (k in c(1:fact_n)) {
        d2q[k,] <- d2[Permut[k,]]
      }
    Bd <- cbind(d1B, d2q)
    Bdq <- Bd[,indUj] # dimensions with respect to permutations
    # Commutator
    M_m_n<-matrix(rep(0,(prod(d1)*prod(d2))^2),nrow = prod(d1)*prod(d2))
      for  (kk in c(1:fact_n)) {
        M_m_n<- M_m_n + .matr_Commutator_Kperm(q[kk,],Bdq[kk,])
      }
  }

  if (useSparse==TRUE){
    i1<- Matrix::Matrix(data =c(1:n,(1:n)+n) , nrow = 2, ncol = n, byrow =TRUE,sparse=TRUE)
    fact_n<-factorial(n)  # number of permutations
    B <- Matrix::Matrix(data = rep(c(1:n),fact_n), nrow = fact_n, ncol = n, byrow =TRUE,sparse=TRUE)
    Permut <- arrangements::permutations(1:n)#perm(c(1:n))
    q<-cbind(B, (Permut+n))
    indUj<- as.vector(i1) # reorder q
    q<-q[,indUj] # permutations
    # dimensions
    d1B <- Matrix::Matrix(rep(d1,fact_n), nrow = fact_n, ncol = n, byrow =TRUE,sparse=TRUE)
    d2q<- Matrix::Matrix(rep(0,fact_n*n),nrow = fact_n, ncol = n,sparse=TRUE) # zeros(fact_n,n);
    for (k in c(1:fact_n)) {
      d2q[k,] <- d2[Permut[k,]]
    }
    Bd <- cbind(d1B, d2q)
    Bdq <- Bd[,indUj] # dimensions with respect to permutations
    # Commutator
    M_m_n<-Matrix::Matrix(rep(0,(prod(d1)*prod(d2))^2),nrow = prod(d1)*prod(d2),sparse=TRUE)
    for  (kk in c(1:fact_n)) {
      M_m_n<- M_m_n + .matr_Commutator_Kperm(q[kk,],Bdq[kk,],useSparse=TRUE)
    }
  }

  return(M_m_n)
}



.indx_Commutator_Mixing <- function(x, d1,d2) {
  if (length(x)!= prod(d1)*prod(d2)) (stop("x must have dimension prod(d1)*prod(d2)"))
  # permutations
  # dim(d1)<-dim(d2)
  n <- length(d2)

  i1<- matrix(data =c(1:n,(1:n)+n) , nrow = 2, ncol = n, byrow =TRUE)
  fact_n<-factorial(n)  # number of permutations
  B <- matrix(data = rep(c(1:n),fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
  Permut <- arrangements::permutations(1:n)#perm(c(1:n))
  q<-cbind(B, (Permut+n))
  indUj<- c(i1) # reorder q
  q<-q[,indUj] # permutations
  # dimensions
  d1B <- matrix(rep(d1,fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
  d2q<- matrix(rep(0,fact_n*n),nrow = fact_n, ncol = n) # zeros(fact_n,n);
  for (k in c(1:fact_n)) {
    d2q[k,] <- d2[Permut[k,]]
  }
  Bd <- cbind(d1B, d2q)
  Bdq <- Bd[,indUj] # dimensions with respect to permutations
  # Commutator
  lcx<-0
  for  (kk in c(1:fact_n)) {
    lcx<- lcx + x[.indx_Commutator_Kperm(q[kk,],Bdq[kk,])]
  }
  return(lcx)
}





############################
### 5

#' Symmetrizer Matrix
#'
#' Based on Chacon and Duong (2015) efficient recursive algorithms for functionals based on higher order
#' derivatives.  An option for sparse matrix is provided. By using sparse matrices  far
#' less memory is required and faster computation times are obtained
#'
#'
#'
#' @param d dimension of a vector x
#' @param n  power of the Kronecker product
#' @param useSparse TRUE or FALSE. If TRUE an object of the class
#' "dgCMatrix" is produced.
#' @return A Symmetrizer matrix with order  \eqn{{d^n} \times d^n}. If \code{useSparse=TRUE}
#' an object of the class "dgCMatrix" is produced.
#'
#' @examples
#' a<-c(1,2)
#' b<-c(2,3)
#' c<-kronecker(kronecker(a,a),b)
#' ## The symmetrized version of c is
#' as.vector(SymMatr(2,3)%*%c)
#'
#' @references Chacon, J. E., and Duong, T. (2015). Efficient recursive algorithms
#' for functionals based on higher order derivatives of the multivariate Gaussian density.
#' Statistics and Computing, 25(5), 959-974.
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.Section 1.3.1 Symmetrization, p.14. (1.29)
#'
#'
#' @family Matrices and commutators
#' @export
SymMatr<-function(d,n,useSparse=FALSE){
  if(useSparse==FALSE) {Id=diag(d)
  if (n==0) {return(1)}
  if (n==1) {return(diag(d))}
  }
  if (useSparse==TRUE) {Id<-Matrix::Diagonal(d)
  if (n==0) {return(1)}
  if (n==1) {return(Id)}
  }
  Symm = Id
  T=Symm
  if (useSparse==FALSE) {A <- .matr_Commutator_Kperm(c(2,1),d)}
  if (useSparse==TRUE) {A <- .matr_Commutator_Kperm(c(2,1),d,useSparse=TRUE)}
  for (k in 2:n){
    T= A%*%kronecker(T,Id)%*%A+A
    Symm = kronecker(Symm,Id)%*%T
    if (k < n) { A=kronecker(Id,A)}
  }

  SymmU = Symm/factorial(n);
  return(SymmU)
}


#' Symmetrizing vector
#'
#' Vector symmetrizing a T-product of vectors of the same dimension d.
#' Produces the same results as SymMatr
#'
#'
#' @param x the vector  to be symmetrized of dimension d^n
#' @param d size of the single vectors in the product
#' @param n power of the T-product
#'
#' @return A vector with the symmetrized version of x of dimension d^n
#'
#' @examples
#' a<-c(1,2)
#' b<-c(2,3)
#' c<-kronecker(kronecker(a,a),b)
#' ## The symmetrized version of c is
#' SymIndx(c,2,3)
#'
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.Section 1.3.1 Symmetrization, p.14. (1.29)
#'
#'
#' @family Matrices and commutators
#' @export

SymIndx  <-  function(x,d,n){
  x.sym <- rep(0,d^n)
  if (d>100)  (stop("d must NOT be greater than 100"))
  xp<- .primnum(545)
  xp<- xp[1:d]
  y <- xp
  for (k in 1:(n-1)){ y <- kronecker(xp,y)}

  di.indx <-match(unique(y), y)
  for (k in 1:length(di.indx)) {
    x.sym[rep(y[di.indx[k]],d^n)==y] <-  mean( x[rep(y[di.indx[k]],d^n)==y] )
  }
  return(x.sym)
}


##########################
### 6

#' Elimination Matrix
#'
#' Eliminates the duplicated/q-plicated elements  in a T-vector of multivariate moments
#' and cumulants.
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section  1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#' @param useSparse TRUE or FALSE.
#' @return Elimination matrix of order \eqn{\eta_{d,q} \times d^q= {d+q-1 \choose q}}.
#' If \code{useSparse=TRUE} an object of the class "dgCMatrix" is produced.
#'
#' @examples
#' x<-c(1,2,3)
#' y<-kronecker(kronecker(x,x),x)
#' ## Distinct elements of y
#' z<-as.matrix(EliminMatr(3,3))%*%y
#' ## Restore eliminated elements in z
#' as.vector(QplicMatr(3,3)%*%z)
#'
#' @family Matrices and commutators
#' @export
EliminMatr<-function(d,q,useSparse=FALSE){
  if (d>100)  (stop("d must NOT be greater than 100"))
  x<- .primnum(545)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  int<-match(unique(y), y)
  if (useSparse==TRUE) {I <- Matrix::Diagonal(d^q)} else {I<- diag(d^q)}
  return(I[int,])
}

###############################
### 7

#' Qplication Matrix
#'
#' Restores the duplicated/q-plicated  elements which are eliminated
#' by  EliminMatr in a T-product of vectors of dimension d.
#'
#'
#' Note: since the algorithm of elimination is not unique, q-plication works together
#' with the function  EliminMatr only.
#'
#'
#'
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, p.21, (1.31)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#' @param useSparse TRUE or FALSE.
#' @return Qplication matrix of order \eqn{d^q \times \eta_{d,q}}, see (1.30), p.15.
#' If \code{useSparse=TRUE} an object of the class "dgCMatrix" is produced.
#'
#' @examples
#' x<-c(1,2,3)
#' y<-kronecker(kronecker(x,x),x)
#' ## Distinct elements of y
#' z<-as.matrix(EliminMatr(3,3))%*%y
#' ## Restore eliminated elements in z
#' as.vector(QplicMatr(3,3)%*%z)
#'
#' @family Matrices and commutators
#' @export
QplicMatr<-function(d,q,useSparse=FALSE){
  if (d>100)  (stop("d must NOT be greater than 100"))
  x<- .primnum(545)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  int<-match(unique(y), y)
  im<-match(y, unique(y))
  esz <- length(int)
  if (useSparse==FALSE) {
    De <- matrix(0,d^q,esz)
    I<-diag(d^q)}
  else {
    De <- Matrix::Matrix(0,d^q,esz)
    I<-Matrix::Diagonal(d^q)}

  for (k in 1:esz){
    if (sum((im==k))>1){De[,k]<-apply(I[ ,(im==k)],1,sum)}
    else {De[,k]<-I[,(im==k)]}
  }
  return(De)
}

#' Qplication vector
#'
#' Restores the duplicated/q-plicated  elements which are eliminated
#' by  EliminMatr or EliminIndx in a T-product of vectors of dimension d.
#' It produces the same results as QplicMatr.
#'
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, p.21, (1.31)
#'
#' @param d dimension of the vectors in the T-product
#' @param q  power of the Kronecker product
#'
#' @examples
#' x<-c(1,2,3)
#' y<-kronecker(kronecker(x,x),x)
#' ## Distinct elements of y
#' z<-y[EliminIndx(3,3)]
#' ## Restore eliminated elements in z
#' z[QplicIndx(3,3)]
#'
#' @return A vector (T-vector) with all elements previously eliminated by EliminIndx
#'
#' @family Matrices and commutators
#' @export
QplicIndx <-function(d,q){
  if (d>100)  (stop("d must NOT be greater than 100"))
  xp<- .primnum(545)
  xp<- xp[1:d]
  y <- xp
  for (k in 1:(q-1)){ y <- kronecker(xp,y)}
  ind.q <- match(y, unique(y))

  return(ind.q)

}





#' Univariate moments and cumulants from T-vectors
#'
#' A vector of indexes to select the moments and cumulants of the single components
#' of the random vector X for which a T-vector of moments and cumulants is available
#'
#' @param d dimension of a vector X
#' @param q  power of the Kronecker product
#' @return A vector of indexes
#'
#' @examples
#' ## For a 3-variate skewness and kurtosis vectors estimated from data, extract
#' ## the skewness and kurtosis estimates for each of the single components of the vector
#' alpha<-c(10,5,0)
#' omega<-diag(rep(1,3))
#' X<-rSkewNorm(200, omega, alpha)
#' EVSK<-SampleEVSK(X)
#' ## Get the univariate skewness and kurtosis for X1,X2,X3
#' EVSK$estSkew[UnivMomCum(3,3)]
#' EVSK$estKurt[UnivMomCum(3,4)]
#' @family Matrices and commutators
#' @export
UnivMomCum<-function(d,q){
  if (d>100)  (stop("d must NOT be greater than 100"))
  x<- .primnum(545)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  ind<-match(x^q, y)  # unique(y)
  return(ind)
}

#' Distinct values selection vector
#'
#' Eliminates the duplicated/q-plicated elements in a T-vector of multivariate moments
#' and cumulants. Produces the same results as EliminMatr.
#' Note EliminIndx does not provide the same results as unique()
#'
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section  1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#'
#' @return A vector of indexes of the distinct elements in the T-vector
#'
#' @examples
#' x<-c(1,0,3)
#' y<-kronecker(x,kronecker(x,x))
#' y[EliminIndx(3,3)]
#' ## Not the same results as
#' unique(y)
#'
#' @family Matrices and commutators
#' @export
EliminIndx<-function(d,q){
  if (d>100)  (stop("d must NOT be greater than 100"))
  x<- .primnum(545)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  ind<-match(unique(y), y)  #
  return(ind)
}




.indx_Commutator_Kmn <-  function(m,n){
  pp <- 1:(m*n)
  #p2 <- t(matrix(data=pp,nrow = m, ncol = n))
  p2<- matrix(data=pp,nrow = n, ncol = m,byrow=TRUE)
  return( as.vector(p2))
}


.indx_Commutator_Kperm <- function(perm,dims){
  n <- length(perm)
  if (length(dims)== 1) {
    perm <- PermutationInv(perm) # included by Gy
    d <- dims
    k <- length(perm)  #
    perm0 <- perm    #
    perm01 <- rev(perm0)
    perm02 <-  order(perm01)
    permutation0 <- rev(perm02)

    pik <- function(perm,d) { 1+sum((perm-1)*cumprod(d*rep(1,k))/d)}

    invIk0 <-  order(permutation0)

    Allind <-  unique(arrangements::combinations( rep(1:d, k), k ) )

    ind1 <- NULL
    ind2 <- NULL
    for (jj in c(1:(d^k))) {
      pA <- Allind[jj,]
      pInv <- pA[invIk0] #inverzepermutation 0
      ind2 <- pik(pA,d)
      ind1[ind2] <- pik(pInv,d)

    }
    return(ind1)   #ind2[ind1]cbind(,ind2)
  }
  S.perm <-  sort(perm , index.return = TRUE)
  perm <- S.perm$ix # % inverse ( sort is increasing)
  D <- prod(dims)
  K.indx <- 1:D #prod(dims)
  if (length(perm)== 1) {return( K.indx )}
  U.a <- NULL
  U.a2 <- NULL
  for (i in 1 : (n - 1)) {
    # run loop (n-i) times
    for (j in 1 : (n - i)) {
      # compare elements
      # j <- 3
      if (perm[j] > perm[j + 1]) {
        if ((j-1)==0) {
          csere <- .indx_Commutator_Kmn(dims[j+1],dims[j])
          dj0 <- dims[j+1]*dims[j]
          dj <- D/dj0
          csere.b <- rep((csere-1)*dj+1, each=dj) + rep(0:(dj-1), dj0)
          K.indx <- K.indx[csere.b]
        }  else {U.a <- prod(dims[1:(j-1)])
        if ((j+1) < n) { U.a2 <- prod(dims[(j+2):n])
        csere <- .indx_Commutator_Kmn(dims[j+1],dims[j])
        dj0.a <- dims[j+1]*dims[j]
        D.a <- U.a*dj0.a
        dj.a <- D.a/dj0.a
        # after 1
        csere.a <- rep(csere, dj.a) + rep(dj0.a*(0:(dj.a-1)), each=dj0.a)
        # K.indx <- K.indx[csere.b]
        # before csere
        dj0.b <- D.a  #D/dims[j+1]/dims[j]
        dj.b <- D/dj0.b  #dims[j+1]*dims[j]
        # K.indx <- rep((csere-1)*dj+1, each=dj) + rep(0:(dj-1), (dims[j+1]*dims[j]))
        csere.b <- rep((csere.a-1)*dj.b+1, each=dj.b) + rep(0:(dj.b-1),dj0.b)

        #   K.indx <- rep(csere, dj) + rep(dj0*(0:(dj-1)), each=dj0)
        K.indx <- K.indx[csere.b]
        } else { U.after <- NULL
        csere <- .indx_Commutator_Kmn(dims[j+1],dims[j])
        dj0 <- dims[j+1]*dims[j]
        dj <- D/dj0
        csere.a <-  rep(csere, dj) + rep(dj0*(0:(dj-1)), each=dj0)
        K.indx <- K.indx[csere.a]
        }
        }
        temp.d <- dims[j]
        dims[j] <- dims[j+1]
        dims[j+1] <- temp.d
        temp.p <- perm[j]
        perm[j] <- perm[j + 1]
        perm[j + 1] <- temp.p
      }
    }
  }
  return(K.indx)
}


.matr_Commutator_Moment<-function(el_rm,d,useSparse=FALSE) {
  N<-length(el_rm)
  PTB<-PartitionTypeAll(N)
  loc_type_el <- .Partition_Type_eL_Location(el_rm)
  r <- loc_type_el[1]
  m <- loc_type_el[2]

  part_class<-PTB$Part.class
  S_N_r<-PTB$S_N_r
  S_m_j<-PTB$S_r_j

  sepL<-cumsum(S_N_r)
  sepS_r<-cumsum(S_m_j[[r]])

  if (useSparse==FALSE){
    if (r==1) {perm_Urk1<- 1:N
    L_eL<-rep(0,d^N)
    L_eL<-L_eL+.matr_Commutator_Kperm(perm_Urk1,d)
    return("LeL"=L_eL)
    }
    else {
      if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
      if (.is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
      else {u_ind<-l_ind+S_m_j[[r]][m]-1}
      perm_Urk1<-matrix(0,S_m_j[[r]][m],N)
      sz<- 1
      for (k in l_ind:u_ind){
        perm_Urk1[sz,]<-.Partition_2Perm(part_class[[k]])
        sz<-sz+1
      }
    }
    L_eL<-rep(0,d^N)
    for (ss in 1:dim(perm_Urk1)[1]) {
      L_eL<-L_eL + .matr_Commutator_Kperm(perm_Urk1[ss,],d)
    }
  }

  if (useSparse==TRUE){
    if (r==1) {perm_Urk1<- 1:N
    L_eL<-rep(0,d^N)
    L_eL<-L_eL+.matr_Commutator_Kperm(perm_Urk1,d,useSparse = TRUE)
    return("LeL"=L_eL)
    }
    else {
      if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
      if (.is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
      else {u_ind<-l_ind+S_m_j[[r]][m]-1}
      perm_Urk1<-Matrix::Matrix(0,S_m_j[[r]][m],N,sparse=TRUE)
      sz<- 1
      for (k in l_ind:u_ind){
        perm_Urk1[sz,]<-.Partition_2Perm(part_class[[k]])
        sz<-sz+1
      }
    }
    L_eL<-rep(0,d^N)
    for (ss in 1:dim(perm_Urk1)[1]) {
      L_eL<-L_eL + .matr_Commutator_Kperm(perm_Urk1[ss,],d,useSparse = TRUE)
    }
  }

  if (useSparse==TRUE) {
    L_eL<-as.matrix(L_eL)
    L_eL<-t(L_eL)
    L_eL<- Matrix::Matrix(L_eL,sparse=TRUE)
    return("L_eL"= L_eL)}

  L_eL<-t(L_eL)

  return("L_eL"= L_eL)
}




.indx_Commutator_Moment<-function(x,el_rm,d) {
  N<-length(el_rm)
  PTB<-PartitionTypeAll(N)
  loc_type_el <- .Partition_Type_eL_Location(el_rm)
  r <- loc_type_el[1]
  m <- loc_type_el[2]

  part_class<-PTB$Part.class
  S_N_r<-PTB$S_N_r
  S_m_j<-PTB$S_r_j

  sepL<-cumsum(S_N_r)
  sepS_r<-cumsum(S_m_j[[r]])


  if (r==1) {perm_Urk1<- 1:N
  px<-0
  px<-px+ x[.indx_Commutator_Kperm(perm_Urk1,d)]

  return("px"=px)
  }
  else {
    if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
    if (.is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
    else {u_ind<-l_ind+S_m_j[[r]][m]-1}
    perm_Urk1<-matrix(0,S_m_j[[r]][m],N)
    sz<- 1
    for (k in l_ind:u_ind){
      perm_Urk1[sz,]<-.Partition_2Perm(part_class[[k]])
      sz<-sz+1
    }
  }
  px<-0
  for (ss in 1:dim(perm_Urk1)[1]) {
    uu<-PermutationInv(perm_Urk1[ss,])
    px<-px+ x[.indx_Commutator_Kperm(uu,d)]
  }


  return("px"= px)
}

