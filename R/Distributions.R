## 1. rUnis
## 2. rSkewNorm
## 3. rCFUSN
## 4. rCFUSSD
## 5. MomCumSkewNorm
## 6. MomCumUnif
## 7. MomCumCFUSN
## 8. MomCumZabs
## 9. distr_UniAbs_EVSK_Th_
## 10. EVSKSkewNorm
## 11. distr_Uni_EVSK_Th
## 12. MomCumMVt - multivariate t
## 13. .evenMomt supports MomCumMVt
## 14. EVSKSkewt
## 15. EVSKGenHyp
## 16. MomCumGenHyp

#####################
#########################################################
#########################    #############
##########################################################
#' Random multivariate spherically symmetric distributions

#'
#' Generate random d-vectors from the spherically symmetric uniform distribution on the sphere
#' @param n sample size
#' @param d dimension
#'
#' @return  A random matrix \eqn{n \times d}
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Random generation
#' @export

rUniS <-  function(n,d){
Z <- matrix( stats::rnorm(d*n), nrow = n)
U3 <- t(apply(Z, 1, function(y) y/sqrt(sum(y^2))))
return(U3)
}


########################################################
###########      ###################
#######################################################
#' Random Multivariate Skew Normal
#'
#' Generate random d-vectors  from the multivariate Skew Normal distribution
#' @param n sample size
#' @param omega correlation matrix   with d dimension
#' @param alpha shape parameter vector of dimension d
#' @return  A random matrix \eqn{n \times d}
#'
#' @examples
#' alpha<-c(10,5,0)
#' omega<-diag(3)
#' x<-rSkewNorm(20,omega,alpha)
#'
#' @references Azzalini, A. with the collaboration of Capitanio, A. (2014).
#' The Skew-Normal and Related Families. Cambridge University Press,
#'   IMS Monographs series.
#' @references Gy.H.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Section 5.1.2
#' @family Random generation
#' @export
rSkewNorm  <- function(n,omega,alpha){
  Om<-omega
  alf<-alpha
  d <- length(alf)
  delt<- as.vector((t(Om%*%alf)))*(1/as.vector(sqrt(1+t(alf)%*%Om%*%alf)))
  SI <- rbind(c(1,t(delt)),cbind(delt,Om))
  ieg<- eigen(SI)
  V <- ieg$vectors
  Sig <- V %*% diag(sqrt(ieg$values)) %*% t(V)
  Z <- matrix(stats::rnorm((d+1)*n),nrow=n)
  Z <- (apply(Z, 1, function(y) Sig%*%y))
  X3 <- t(Z[2:(d+1),])*
    matrix(kronecker(rep(1,d),Z[1,]>0)-kronecker(rep(1,d),Z[1,]<0),nrow=n)
  return(X3)
}



###################################################################
#########################(rCFUSN) ################
#######################################################################

#' Random multivariate CFUSN
#'
#' Generate random d-vectors  from the multivariate
#' Canonical Fundamental Skew-Normal (CFUSN) distribution
#' @param n The number of variates to be generated
#' @param Delta Correlation matrix, the skewness matrix Delta
#' @return  A random matrix \eqn{n \times d}
#' @examples
#' d <- 2; p <- 3
#' Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
#' ieg<- eigen(diag(p)+t(Lamd)%*%Lamd)
#' V <- ieg$vectors
#' Delta <-Lamd %*% V %*% diag(1/sqrt(ieg$values)) %*% t(V)
#' x<-rCFUSN(20,Delta)
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 (5.5) p.247
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Random generation
#' @export

rCFUSN  <- function(n,Delta){
  #
  d <- dim(Delta)[1]
  p <- dim(Delta)[2]
  iegD<- eigen(diag(d)-Delta%*%t(Delta))
  VD <- iegD$vectors

  Delta2 <-VD %*% diag(sqrt(iegD$values)) %*% t(VD)
  Z1 <- abs(matrix(stats::rnorm(p*n), nrow = n))
  Z2 <- abs(matrix(stats::rnorm(d*n), nrow = n))

  X2 <- apply(Z2, 1, function(y) Delta2%*%(y))
  #Delta2%*%Z2[1,]
  X <- Delta%*%t(Z1)+X2
  return(t(X))
}


#####################################################################
#################### Canonical Fundamental Skew-Spherical Distribution CFUSSD
###############################################################


################# #################
#' Random multivariate CFUSSD
#'
#' Generate random d-vectors  from the multivariate Canonical Fundamental
#' Skew-Spherical distribution (CFUSSD)  with Gamma generator
#' @param n sample size
#' @param d  dimension
#' @param p  dimension of the first term of (5.5)
#' @param a  shape  parameter of the Gamma generator
#' @param b  scale parameter of the Gamma generator
#' @param Delta skewness matrix
#'
#' @return  A  matrix of \eqn{n \times d} random numbers
#'
#' @examples
#' n <- 10^3; d <- 2; p <- 3 ; a <- 1; b <- 1
#' Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
#' ieg<- eigen(diag(p)+t(Lamd)%*%Lamd)
#' V <- ieg$vectors
#' Delta <-Lamd %*% V %*% diag(1/sqrt(ieg$values)) %*% t(V)
#' rCFUSSD(20,d,p,1,1,Delta)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, (5.36) p. 266, (see p.247 for Delta)
#' @family Random generation
#' @export
rCFUSSD  <- function(n,d,p,a,b,Delta){
  #
  ieg2 <- eigen(diag(d)-Delta%*%t(Delta))
  V2 <- ieg2$vectors
  Delta2 <-  V2 %*% diag(sqrt(ieg2$values)) %*% t(V2)
  R <- stats::rgamma(n, shape= a, scale = b ) # a > 0 scale = 1/rate
  p1 <-  stats::rbeta(n,p/2,d/2,n)

  U2 <-  rUniS(n,d)
  mU1 <- abs(rUniS(n,p));
  Z1 <-  mU1*matrix(kronecker(rep(1,p),R*p1),ncol = p)
  Z2 <-  U2*matrix(kronecker(rep(1,d),R*(1-p1)),ncol = d)
  X <- t(Delta%*%t(Z1)+Delta2%*%t(Z2))
  return(X)
}



.distr_Uni_EVSK_Th <- function(d,nCum = TRUE){
  eL <- c(0,2,0,0)
  # loc_type_el <- Partition_Type_eL_Location(eL)
  #L2 <- .Commutator_Moment_eL(eL,d)
  Idv <- as.vector(diag(d))
  EU1 <- rep(0,d)
  EU2 <- 1/d*Idv
  varU<- matrix(1/d*Idv,nrow = d)
  EU3 <- rep(0,d^3)
  x<- kronecker(Idv,Idv)
  EU4 <- .indx_Commutator_Moment_t(x,eL,d)*1/d/(d+2)
  EU.k <- list(EU1,EU2,EU3,EU4)

  UEVSK <- EU.k
  names(UEVSK) <- c("EU1","EU2","EU3","EU4")
  if (nCum == TRUE){
    cumU.k <- .conv_Mom2CumMulti(EU.k)
    kurtU <- as.vector(cumU.k[[4]]*d^2)
    skewU <- rep(0,d^3)
    Skew.tot<-c(0)
    Kurt.tot<-sum(kurtU^2)
    UEVSK <- list(EU1,varU,skewU,Skew.tot,kurtU,Kurt.tot)
    names(UEVSK) <- c("EU1","varU","Skew.U","Skew.tot","Kurt.U", "Kurt.tot")
  }
  return(UEVSK)
}


#########################################################################
#' Moments and cumulants Uniform Distribution on the Sphere
#'
#' By default, only moments are provided
#'
#' @param r highest order of moments and cumulants
#' @param nCum if it is TRUE then cumulants are calculated
#' @param d dimension
#' @return The list of moments and cumulants in vector form
#' @examples
#' # The first four moments for d=3
#' MomCumUniS(4,3,nCum=0)
#' # The first four moments and cumulants for d=3
#' MomCumUniS(4,3,nCum=4)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 Proposition 5.3 p.297
#' @family Moments and cumulants
#' @export

MomCumUniS <- function(r,d,nCum = FALSE)
{
  m0 <- floor(r/2)
  # 1d moments
  EU.k0 <- c(rep(0,r))
  EU.k0[seq(2,2*m0,2)] <- cumprod(seq(1,2*m0,2))/cumprod(seq(d,d+2*(m0-1),2))
  Idv <- as.vector(diag(d))
  per <- cumprod(seq(1,2*m0,2))
  Idv.u <- Idv
  EUM.k <- NULL
  EUM.k[[1]]  <- rep(0,d)
  EUM.k[[2]] <- 1/d*Idv
  EUM.k[[3]] <- c(rep(0,d^3))
  varU<- matrix(1/d*Idv,nrow = d)

  if (m0 >=2) {
  for (k in 2:m0) {
    # tic
    eL <- c(0,k,c(rep(0,2*k-2)))
    Idv.u <- kronecker(Idv.u,Idv)
    EUM.k[[2*k]] <- as.vector(EU.k0[2*k] *.indx_Commutator_Moment_t(Idv.u,eL,d)/per[k])
    EUM.k[[2*k+1]] <-c(rep(0,d^(2*k+1)))
    #   print( toc)
  }
  }
  if ( r%%2 == 0 ) EUM.k <- EUM.k[-(2*m0+1)]
  if (nCum == TRUE) {
    CumU <- .conv_Mom2CumMulti(EUM.k )
    MomCumU <- list(EUM.k,CumU )
    names(MomCumU) <- c("EUM" , "CumU")
    return(MomCumU)
  }

  nev <- NULL
  for (k in 1:r) {
    nev[[k]] <- paste("EU.",as.character(k), sep='')
  }
  names(EUM.k) <- nev
  return(EUM.k)
}



.distr_UniAbs_EVSK_Th <- function(d,nCum=TRUE ){

  Id <- diag(d)
EabsU1 <- sqrt(1/pi)/.Gkd(1,d)*rep(1,d)
EabsU2 <- as.vector(Id)/d+(rep(1,d^2)-as.vector(Id))/pi/.Gkd(2,d)
###########
e3 <- rep(0,d^3)
for (k in c(1:d)) {
  e3 <- e3 + .kron3(Id[,k])
}
e3 <- matrix(e3,ncol = 1)
e21 <- rep(0,d^3)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e21 <- e21 +  kronecker(.kron2(Id[,k]),Id[,j])
  }
}

e21 <- SymIndx(e21,d,3)*3
e111 <- rep(0,d^3)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          e111 <- e111 +  kronecker(Id[,k],kronecker(Id[,j],Id[,m]))
        }
      }
    }
  }
}
e111<-matrix(e111,ncol=1) # 6*matr_Symmetry(d,3) %*% e111
EabsU3 <-as.vector(( e3 + e21/2 + e111/pi)/sqrt(pi)/.Gkd(3,d))
#########################################################################
e4 <- rep(0,d^4)
for (k in c(1:d)) {
  e4 <- e4 + .kron4(Id[,k])
}
e4 <-matrix(e4,ncol = 1)
##############################################################x
e31 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e31 <- e31 +  kronecker(.kron3(Id[,k]),Id[,j])
  }
}
e31 <- SymIndx(e31,d,4)*4
###########################################
e211 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          e211 <- e211 +  kronecker(.kron2(Id[,k]),kronecker(Id[,j],Id[,m]))
        }
      }
    }
  }
}

e211 <- SymIndx(e211,d,4)*6
##################################
e22 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e22 <- e22 + kronecker(.kron2(Id[,k]),.kron2(Id[,j]))
  }
}
e22 <- matrix(e22,ncol = 1)
###########################################
e1111 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          for (n in c(1:d)) {
            if (n !=k & n!= j & n!=m) {
              e1111 <- e1111 + kronecker(kronecker(Id[,k],Id[,j]),kronecker(Id[,m],Id[,n]))}
          }

        }
      }
    }
  }
}

e1111 <- SymIndx(e1111,d,4)
###########################################################
e.vect4 <-( e31/pi +  e22/4 + e211/sqrt(pi^3)  +  e1111/pi^2)/.Gkd(4,d)
EabsU4 <- as.vector(3*e4/d/(d+2) +  e.vect4)
Eabs <- list(EabsU1,EabsU2,EabsU3,EabsU4 )

if (nCum == TRUE){
   cumU.k <- .conv_Mom2CumMulti(Eabs)
# Eabs[[1]] <- c(rep(0,d))
varU <- matrix(cumU.k[[2]],nrow = d)
# varU.inv <-  solve(varU)  # %*%varU
U.V <- eigen(varU)$vectors
U.L <- eigen(varU)$values
sqrtU.inv <- U.V %*% diag(1/sqrt(U.L)) %*% t(U.V)
szor.skew <- .kron3(sqrtU.inv )
skew.U <- as.vector(szor.skew%*%cumU.k[[3]])
kurt.U <- as.vector(.kron4(sqrtU.inv) %*% cumU.k[[4]])
    EabsT <- list(EabsU1,varU,Eabs,cumU.k,skew.U,kurt.U)
  names(EabsT) <- c("EU1","varU","EU.k","cumU.k","skew.U","kurt.U")
  return(EabsT)
 }
names(Eabs) <- c("EabsU1","EabsU2","EabsU3","EabsU4" )
return(Eabs)
}




.distr_Zabs_MomCum_Th <- function(r,nCum=FALSE){
   #
   n<-r
   n0 <- floor((n-2)/2) +(n/2-n%/%2>0)
   mu <- NULL
   mu[1] <- sqrt(2/pi)
   mu[2] <- 1
   double.fact1 <- NULL
   double.fact2 <- NULL
   double.fact1[1] <- 1
   double.fact2[1] <- 1
   #sz <- 1
   if (n0>0){
     for (sz in 1:n0) {
       double.fact1[sz+1] <- (2*sz)*double.fact1[sz]
       mu[2*sz+1] <- mu[1]*double.fact1[sz+1]
       double.fact2[sz+1] <- (2*sz+1)*double.fact2[sz]
       mu[2*sz+2] <- double.fact2[sz+1]
     }
   }
   if (n/2-n%/%2 >0) mu <- mu[-(2*n0+2)]
   if (nCum==1){
     cum <- .conv_Mom2Cum(mu )
     MomCumZ <- list(mu,cum )
     names(MomCumZ) <- c("MuZ" , "CumZ")
     return(MomCumZ)}
   return(mu)
 }

##########################

####################################################
#######
########################################################
###########       SkewNormEVSK_Th ###################
#######################################################
 #' EVSK multivariate Skew Normal
 #'
 #' Computes the theoretical values of the mean vector,  covariance, skewness vector,
 #' total skenwness, (excess) kurtosis vector and total
 #' kurtosis for the multivariate Skew Normal distribution
 #' @param omega A \eqn{d \times d} correlation matrix
 #' @param alpha shape parameter d-vector
 #' @return  A list of theoretical values for the mean vector,  covariance, skewness vector,
 #' total skenwness, kurtosis vector and total kurtosis
 #'
 #' @examples
 #' alpha<-c(10,5,0)
 #' omega<-diag(3)
 #' EVSKSkewNorm(omega,alpha)
 #' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
 #' Springer 2021 (5.5) p.247
 #' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
 #' skewness and kurtosis. Sankhya A, 83(2), 607-644.
 #' @family Moments and cumulants
 #' @export

EVSKSkewNorm <- function(omega,alpha){

  oszt <- sqrt(1+t(alpha)%*%omega%*%alpha)
  dim(oszt) <- NULL
  delt <- omega%*%alpha/oszt
  EX <- as.vector(sqrt(2/pi)*delt)
  VarX <- omega-2/pi*delt%*%t(delt)
  ieg <- eigen(VarX)
  V <- ieg$vectors
  VarX.inv <-  V %*% diag(1/sqrt(ieg$values)) %*% t(V)
  VarXid <- VarX.inv %*% delt
  kr3 <- .kron3(VarXid)
  SkewX <-as.vector((2*sqrt(2/pi)^3-sqrt(2/pi))*kr3)
  SkewX.tot <- sum(SkewX*SkewX)
  kz4  <-  -6*4/(pi^2)+4*2/pi
  KurtX <- as.vector(kz4*kronecker(kr3,VarXid))
  KurtX.tot <- sum(KurtX*KurtX)

  EVSKTheo <- list(EX,VarX, SkewX, SkewX.tot,KurtX, KurtX.tot )
  names(EVSKTheo) <- c("EX","VarX", "SkewX", "SkewX.tot","KurtX", "KurtX.tot")
  return(EVSKTheo)
}

########################################################
###########       MomCumSkewNorm ###################
#######################################################
#' Moments and cumulants d-variate Skew Normal
#'
#'
#' Computes the theoretical values of moments and cumulants up to the r-th order.
#' Warning: if nMu = TRUE it can be very slow
#' @param r the highest moment and cumulant order
#' @param omega A \eqn{d \times d} correlation matrix
#' @param alpha shape parameter d-vector
#' @param nMu if it is TRUE then moments are calculated as well
#' @return  A list of theoretical moments and cumulants
#'
#' @examples
#' alpha<-c(10,5,0)
#' omega<-diag(3)
#' MomCumSkewNorm(r=4,omega,alpha)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 (5.5) p.247, Lemma 5.1 p. 246
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Moments and cumulants
#' @export

MomCumSkewNorm <- function(r = 4,omega,alpha,nMu = FALSE ){
  #
  # Om correlation matrix
  Om<-omega
  alf<-alpha
  oszt <- sqrt(1+t(alf)%*%Om%*%alf)
  dim(oszt) <- NULL
  delt <- Om%*%alf/oszt
  EX <- sqrt(2/pi)*delt
  VarX <- Om-2/pi*delt%*%t(delt)
  deltk <- NULL
  deltk[[1]] <- delt
  deltk[[2]] <- kronecker(delt,delt)
  MC.absZ <- .distr_Zabs_MomCum_Th(r,nCum = 1)
  cum <- NULL
  cum[[1]] <- as.vector(EX)
  cum[[2]] <- as.vector(VarX)
  for (k in 3:r) {
    deltk[[k]] <- kronecker(delt,deltk[[k-1]])
    cum[[k]] <- as.vector(MC.absZ$CumZ[k]* deltk[[k]])
  }
 if (nMu == TRUE){
  Mu <- .conv_Cum2MomMulti(cum)
  EVSKTheo <- list(Mu,cum )
  names(EVSKTheo) <- c("MuX","CumX")
  return(EVSKTheo)
 }
  EVSKTheo <- list(cum )
  names(EVSKTheo) <- c("CumX")
  return(EVSKTheo)

}



.distr_ZabsM_MomCum_Th  <- function(r,d,nCum=FALSE){
  #
  MC.Z <- .distr_Zabs_MomCum_Th(r,nCum = 1)
  Id <- diag(d)
  i_kron_k <- NULL
  i_kron_k[[1]] <- rep(1,d)
  e_kron_k <- NULL
  e_kron_k[[1]] <-Id
  CumMZ <- NULL
  for (k in 1:r) {
    i_kron_k[[k+1]] <- rep(0,d^(k+1))
    #e_kron_k[[k+1]] <- NULL
    e_kron_k[[k+1]] <- matrix(rep(0,d^(k+2)),ncol=d)
    for (j in 1:d) {
      e_kron_k[[k+1]][,j] <- kronecker(e_kron_k[[k]][,j],Id[,j])
      i_kron_k[[k+1]] <- i_kron_k[[k+1]]+e_kron_k[[k+1]][,j]
    }
    CumMZ[[k]] <-  MC.Z$CumZ[k]*i_kron_k[[k]]
  }
  MuMZ <- .conv_Cum2MomMulti(CumMZ)
  if (nCum == TRUE){
    MomCumMZ <- list(MuMZ,CumMZ )
    names(MomCumMZ) <- c("MuMZ","CumMZ")
    return(MomCumMZ)
  }
  return(MuMZ)
}


#############################################
#' Moments and cumulants CFUSN
#'
#' Provides the theoretical cumulants of the multivariate Canonical Fundamental
#' Skew Normal distribution
#' @param r The highest cumulant order
#' @param d The multivariate dimension and number of rows of the skewness matrix Delta
#' @param p The number of cols of the skewness matrix Delta
#' @param Delta The skewness matrix
#' @param nMu If set to TRUE, the list of the first r d-variate moments is provided
#' @return  The list of theoretical  cumulants in vector form
#' @examples
#' r <- 4; d <- 2; p <- 3
#' Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
#' ieg<- eigen(diag(p)+t(Lamd)%*%Lamd)
#' V <- ieg$vectors
#' Delta <-Lamd %*% V %*% diag(1/sqrt(ieg$values)) %*% t(V)
#' MomCum <- MomCumCFUSN(r,d,p,Delta)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Lemma 5.3 p.251
#' @family Moments and cumulants
#' @export

MomCumCFUSN  <- function(r,d,p,Delta,nMu = FALSE){
  #
  EX <- sqrt(2/pi)*Delta%*%rep(1,p)
  VarX <- diag(d)-2/pi*Delta%*%t(Delta)
  MC.Z <- .distr_Zabs_MomCum_Th(r,nCum = 1)
  Ip <- diag(p)
  i_kron_k <- NULL
  i_kron_k[[1]] <- rep(1,p)
  e_kron_k <- NULL
  e_kron_k[[1]] <- Ip
  De_kron_k <- NULL
  De_kron_k[[1]] <- Delta
  CumMX <- NULL
  CumMX[[1]] <- as.vector(EX)
  if (r==1) {return(CumMX)}
  for (k in 1:(r-1)) {
    i_kron_k[[k+1]] <- rep(0,p^(k+1))
    # e_kron_k[[k+1]] <- NULL
    e_kron_k[[k+1]] <- matrix(rep(0,p^(k+2)),ncol=p)
    for (j in 1:p) {
      e_kron_k[[k+1]][,j] <- kronecker(e_kron_k[[k]][,j],Ip[,j])
      i_kron_k[[k+1]] <- i_kron_k[[k+1]]+e_kron_k[[k+1]][,j]
    }
    De_kron_k[[k+1]] <- kronecker(De_kron_k[[k]],Delta)
    CumMX[[k+1]] <-  as.vector(MC.Z$CumZ[k]*De_kron_k[[k+1]]%*%(i_kron_k[[k+1]]))
  }
  if (nMu == TRUE){
    MuMX <- .conv_Cum2MomMulti(CumMX)
    MomCumMX <- list(MuMX,CumMX )
    names(MomCumMX) <- c("MuMX","CumMX")
    return(MomCumMX)
  }
  return(CumMX)
}


#' Moments and cumulants Multivariate t-Student distribution
#'
#' The t- distribution is defined as
#' \deqn{\mathbf{X} = \sqrt{\frac{p}{S^2}} \mathbf{Z}}
#' where \eqn{\mathbf{Z}} is a multivariate standard-normal random vector
#' and  \eqn{S^2} is a \eqn{\chi^2_p}.
#' random variable independent of \eqn{\mathbf{Z}}.
#' @param p degrees of freedom
#' @param d dimension
#' @param r highest order of moments and cumulants
#' @param nCum if it is TRUE then cumulants are calculated

#' @return The list of moments (or cumulants) in vector form
#' @examples
#' # The first four moments for trivariate t with 10 d.f.
#' MomCumMVt(p=10,d=3,r=4,nCum=FALSE)
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 Proposition XXXXXX
#' @family Moments and cumulants
#' @export

MomCumMVt <- function(p, d, r,nCum=FALSE) {
  # p: degrees of freedom
  # d: vector dimension
  # r: highest moment order (up to which we compute)
  if (d<2) {stop("Error: 'd' should be greater than 1")}
  moments_list <- list()

  for (k in 1:r) {
    if (k %% 2 == 0) {
      # Even-order moment (use momMt)
      m <- k / 2
      moments_list[[paste0("mu", k)]] <- .evenMomt(p, d, m)
    } else {
      # Odd-order moment is zero vector of length d^k
      zero_vec <- rep(0, d^k)
      moments_list[[paste0("mu", k)]] <- zero_vec
    }
  }

  if (nCum == TRUE) {moments_list <- Mom2Cum(moments_list,Type="Multivariate")}

  return(moments_list)
}


# Computes even moments for the multivariate t distribution
# with p degrees of freedom
.evenMomt <- function(p,d,m){
  # p - degreees of freedom
  # d - vector dimension
  # 2*m even order moment
  if (p/2 <= m ) {stop("Error: 'p' should be greater than 2m")}
  C1 <- p^m/(2^m * gamma(p/2)/gamma(p/2-m))
  const <- C1 * .double_factorial(2*m-1)
  Im <- .kron_power(c(diag(d)),m)
  demom <- const * SymIndx(Im,d,(2*m))
  return(demom)
}


#' EVSK multivariate Skew-t
#'
#' Computes the theoretical values of the mean,  variance,
#' skewness and (excess) kurtosis vectors for the d-variate Skew-t distribution \eqn{St_d(\xi, \boldsymbol{\Omega},
#' \boldsymbol{\alpha},m)}
#' defined as
#'  \deqn{Y = \xi + \sqrt{\frac{m}{S^2}} \mathbf{X}}
#' where \eqn{\mathbf{X}} is a multivariate skew-normal random variable
#' \eqn{SN_d(0, \boldsymbol{\Omega} , \boldsymbol{\alpha})} and \eqn{S^2} is a \eqn{\chi^2_m}
#' random variable independent of \eqn{\mathbf{X}}.
#' @param xi A mean vector
#' @param omega A \eqn{d \times d} correlation matrix
#' @param alpha shape parameter d-vector
#' @param m degrees of freedom
#' @return  A list of theoretical values for the mean,   variance, skewness and
#' kurtosis vectors
#'
#' @examples
#' xi <- c(0,0,0) #
#' alpha <- c(10,5,0) #
#' omega <- diag(3) #
#' m <- 10 #
#'
#' EVSKSkewt(xi,omega,alpha,m)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021  p.277
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Moments and cumulants
#' @export
EVSKSkewt <- function(xi,omega,alpha,m){
  nxi <- xi #
  nomega <- omega #
  nalpha <- alpha #
    ndf <- m #
  d <- length(alpha)
  ndelta <- c(nalpha%*%nomega)/c((1+nalpha%*%nomega%*%nalpha)^(1/2))
  #########  mean and variance
  # mean vector
  muSt <- nxi + sqrt(ndf/pi) * ndelta * .Gkd(-1,ndf)
  # variance vector
  K2St <- (ndf/(ndf-2))*c(nomega)- (ndf/pi) * .Gkd(-1,ndf)^2 * c(kronecker(ndelta,ndelta))
  # variance matrix
  K2m <- matrix(K2St,d,d)
  # SVD of K2m
  svdx<-svd(K2m)
  K2_m12<-svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  ### Skewness terms
  C1 <- sqrt(2/pi) * (ndf/2)^(3/2) * .Gkd(-1,ndf)
  C2 <- (4/pi) * .Gkd(-1,ndf)^2 - 2 / (ndf-3)
  C3 <-  2/((ndf-3)*(ndf-2))
  D1 <- c(.kron_power(ndelta,3) )
  D2 <- c(kronecker(c(nomega),ndelta))
  D2 <- D2 + D2[CommutatorIndx(Type="Kperm",perm=c(1,3,2),dims=d)] + D2[CommutatorIndx(Type="Kperm",perm=c(3,1,2),dims=d)]
  cum3St <-  C1 * (C2 * D1 + C3 * D2 )
  skewSt  <- c(.kron_power(K2_m12,3)%*%cum3St)
  ### kurtosis terms
  KC1 <- (4/pi) * (ndf/2)^2 * .Gkd(-1,ndf)^2 *(4/(ndf-3) -(6/pi) * .Gkd(-1,ndf)^2 )
  KC2 <- (ndf/2)^2 * 8 / ( (ndf-4) * (ndf-2 )^2 )
  KC3 <- (8/pi) * (ndf/2)^2 * .Gkd(-1,ndf)^2 / ( (ndf-3) * (ndf-2 ) )
  KD1 <- .kron_power(ndelta,4)
  KD2 <- .kron_power(c(nomega),2)
  KD2 <- KD2 + KD2[CommutatorIndx(Type="Kperm",perm=c(1,3,2,4),dims=d)] + KD2[CommutatorIndx(Type="Kperm",perm=c(1,4,2,3),dims=d)]
  KD3 <- c(kronecker(c(nomega),kronecker(ndelta,ndelta)))
  KD3 <- KD3 + KD3[CommutatorIndx(Type="Kperm",perm=c(1,3,2,4),dims=d)] +
    KD3[CommutatorIndx(Type="Kperm",perm=c(3,1,2,4),dims=d)] +
    KD3[CommutatorIndx(Type="Kperm",perm=c(3,1,4,2),dims=d)] +
    KD3[CommutatorIndx(Type="Kperm",perm=c(1,3,4,2),dims=d)] +
    KD3[CommutatorIndx(Type="Kperm",perm=c(3,4,1,2),dims=d)]
  cum4St <- KC1 * KD1 + KC2 * KD2 - KC3 * KD3
  kurtSt <- c(.kron_power(K2_m12,4)%*%cum4St)
  EVSKSt <- list(muSt,K2m, skewSt, kurtSt )
  names(EVSKSt) <- c("mu","variance", "skewness", "kurtosis")
  return(EVSKSt)
}


#' EVSK multivariate Generalized hyperbolic
#'
#' Computes the theoretical values of the mean,  variance,
#' skewness and (excess) kurtosis vectors for the d-variate Generalized
#' Hyperbolic distribution \eqn{\mathcal{GH}\left( \lambda
#' ,\chi ,\psi ,\boldsymbol{\mu },\boldsymbol{\Sigma },\boldsymbol{\gamma }%
#' \right)}
#' defined as
#'  \deqn{\mathbf{X}=\boldsymbol{\mu }+V\boldsymbol{\gamma }+\sqrt{V}\boldsymbol{%
#' \Sigma }^{1/2}\mathbf{Z}}
#' where \eqn{\mathbf{Z}\in \mathcal{N}\left( 0,\mathbf{I}_{d}\right)},
#' \eqn{ V \geq 0}, is independent of \eqn{\mathbf{Z}}, is a non-negative,
#' scalar-valued variate, which is \emph{Generalized Inverse Gaussian} (scalar
#'  valued GIG), \eqn{V\in GIG\left( \lambda ,\chi ,\psi \right)}.
#'
#' @param lambda scalar valued
#' @param chi scalar valued
#' @param psi scalar valued
#' @param mu a vector of dimension d
#' @param sigma a dxd covariance matrix
#' @param gamma a scalar value
#'
#' @return  A list of theoretical values for the mean,   variance, skewness and
#' kurtosis vectors
#'
#' @examples
#' lambda <- 1
#' chi <- 2
#' psi <- 2
#' mu <- rep(0,2)
#' sigma <- diag(2)
#' gamma <-  c(0.2,0.5)
#' EVSKGenHyp(lambda, chi, psi, mu, sigma, gamma)
#' @references  A.J. McNeil, R. Frey, and P. Embrechts.
#' Quantitative risk management: concepts,
#' techniques and tools-revised edition. Princeton university press, 2015.
#' @family Moments and cumulants
#' @export
EVSKGenHyp <- function(lambda, chi, psi, mu, sigma, gamma){

  # --- Coerce/check input shapes  ---
  mu    <- as.numeric(mu)
  gamma <- as.numeric(gamma)

  # ensure sigma is a numeric matrix (square)
  sigma <- as.matrix(sigma)
  if (!is.numeric(sigma) || nrow(sigma) != ncol(sigma)) {
    stop("`sigma` must be a numeric square matrix")
  }

  d <- length(gamma)
  if (length(mu) != d) {
    stop("length(mu) must equal length(gamma)")
  }
  if (nrow(sigma) != d) {
    stop("sigma must be a d x d matrix where d = length(gamma)")
  }

  # Function for E[V^r]
  E_Vr <- function(r, lambda, chi, psi) {
    (chi / psi)^(r / 2) * (besselK(sqrt(chi * psi), lambda + r) / besselK(sqrt(chi * psi), lambda))
  }

  mom_V <- c(E_Vr(1, lambda, chi, psi),
             E_Vr(2, lambda, chi, psi),
             E_Vr(3, lambda, chi, psi),
             E_Vr(4, lambda, chi, psi))

  cum_V <- Mom2Cum(mom_V,Type="Univariate")

  ## - mean GH
  mean_GH <- mu + cum_V[1]*gamma

  ## - Variance GH
  sigma_GH <- cum_V[1] * sigma + cum_V[2] * outer(gamma,gamma)

  ## Third cumulant
  Cum3_GH <- c(SymMatr(d,3) %*% (3 * cum_V[2] * kronecker(gamma,c(sigma)) +
                                   cum_V[3] * kronecker(gamma,kronecker(gamma,gamma))))

  ## Fourth cumulant

  Cum4_GH <- c(SymMatr(d,4) %*% (3 * cum_V[2] * kronecker(c(sigma),c(sigma)) +
                                   6 * cum_V[3] * kronecker(c(sigma),kronecker(gamma,gamma)) +
                                   cum_V[4] * kronecker(gamma,kronecker(gamma,kronecker(gamma,gamma)))   ))

  # SVD of Sigma_GH
  svdx<-svd(sigma_GH)
  S_m12 <- svdx$u %*% diag(1/sqrt(svdx$d)) %*% t(svdx$u)

  ### Skewness
  SkewGH  <- c(.kron_power(S_m12,3) %*% Cum3_GH)

  ### kurtosis terms
  KurtGH <- c(.kron_power(S_m12,4) %*% Cum4_GH)

  EVSKGH <- list(mean_GH,sigma_GH, SkewGH, KurtGH )
  names(EVSKGH) <- c("mean","variance", "skewness", "kurtosis")
  return(EVSKGH)
}


#' Moments and cumulants of the multivariate
#' Generalized Hyperbolic distribution
#'
#' Computes cumulants and moments up to order r=6 of the d-variate Generalized
#' Hyperbolic distribution \eqn{\mathcal{GH}\left( \lambda
#' ,\chi ,\psi ,\boldsymbol{\mu },\boldsymbol{\Sigma },\boldsymbol{\gamma }%
#' \right)}
#' defined as
#'  \deqn{\mathbf{X}=\boldsymbol{\mu }+V\boldsymbol{\gamma }+\sqrt{V}\boldsymbol{%
#' \Sigma }^{1/2}\mathbf{Z}}
#' where \eqn{\mathbf{Z}\in \mathcal{N}\left( 0,\mathbf{I}_{d}\right)},
#' \eqn{ V \geq 0}, is independent of \eqn{\mathbf{Z}}, is a non-negative,
#' scalar-valued variate, which is \emph{Generalized Inverse Gaussian} (scalar
#'  valued GIG), \eqn{V\in GIG\left( \lambda ,\chi ,\psi \right)}.
#'
#' @param r highest order of moments and cumulants
#' @param lambda scalar valued
#' @param chi scalar valued
#' @param psi scalar valued
#' @param mu a vector of dimension d
#' @param sigma a dxd covariance matrix
#' @param gamma a scalar value
#' @param nMu if it is TRUE then moments are calculated
#'
#' @return The list of moments (or cumulants) in vector form
#' @examples
#' lambda <- 1
#' chi <- 2
#' psi <- 2
#' mu <- rep(0,2)
#' sigma <- diag(2)
#' gamma <-  c(0.2,0.5)
#' MomCumGenHyp(r=4,lambda, chi, psi, mu, sigma, gamma)
#' @family Moments and cumulants
#' @export
MomCumGenHyp <- function(r = 4, lambda, chi, psi, mu, sigma, gamma, nMu = FALSE) {

  # --- Coerce/check input shapes  ---
  mu    <- as.numeric(mu)
  gamma <- as.numeric(gamma)

  # ensure sigma is a numeric matrix (square)
  sigma <- as.matrix(sigma)
  if (!is.numeric(sigma) || nrow(sigma) != ncol(sigma)) {
    stop("`sigma` must be a numeric square matrix")
  }

  d <- length(gamma)
  if (d < 2) stop("Error: 'd' should be greater than 1")

  if (length(mu) != d) {
    stop("length(mu) must equal length(gamma)")
  }
  if (nrow(sigma) != d) {
    stop("sigma must be a d x d matrix where d = length(gamma)")
  }

  if (r > 6) {
    warning("Only moments/cumulants up to order 6 are provided.")
    r <- 6
  }

  # Function for E[V^r]
  E_Vr <- function(r, lambda, chi, psi) {
    (chi / psi)^(r / 2) * (besselK(sqrt(chi * psi), lambda + r) / besselK(sqrt(chi * psi), lambda))
  }

  # Precompute up to 6th order only
  mom_V <- sapply(1:6, function(k) E_Vr(k, lambda, chi, psi))
  cum_V <- Mom2Cum(mom_V, Type = "Univariate")

  ## - mean GH
  mean_GH <- mu + cum_V[1] * gamma

  ## - Second cumulant - vec(Variance GH)
  Cum2_GH <- c(cum_V[1] * sigma + cum_V[2] * outer(gamma, gamma))

  ## Third cumulant
  Cum3_GH <- c(SymMatr(d, 3) %*% (
    3 * cum_V[2] * kronecker(gamma, c(sigma)) +
      cum_V[3] * kronecker(gamma, kronecker(gamma, gamma))
  ))

  ## Fourth cumulant
  Cum4_GH <- c(SymMatr(d, 4) %*% (
    3 * cum_V[2] * kronecker(c(sigma), c(sigma)) +
      6 * cum_V[3] * kronecker(c(sigma), kronecker(gamma, gamma)) +
      cum_V[4] * kronecker(gamma, kronecker(gamma, kronecker(gamma, gamma)))
  ))

  cumulant_list <- list(mean_GH, Cum2_GH, Cum3_GH, Cum4_GH)

  if (r >= 5) {
    sigam <- sigma %*% gamma
    Cum5_GH <- c(SymMatr(d, 5) %*% (
      15 * cum_V[3] * kronecker(.kron_power(c(sigma), 2), sigam) +
        10 * cum_V[4] * kronecker(c(sigma), .kron_power(sigam, 3)) +
        cum_V[5] * .kron_power(sigam, 5)
    ))
    cumulant_list[[5]] <- Cum5_GH
  }

  if (r >= 6) {
    sigam <- sigma %*% gamma
    Cum6_GH <- c(SymMatr(d, 6) %*% (
      15 * cum_V[3] * .kron_power(c(sigma), 3) +
        45 * cum_V[4] * kronecker(.kron_power(c(sigma), 2), .kron_power(sigam, 2)) +
        15 * cum_V[5] * kronecker(c(sigma), .kron_power(sigam, 4)) +
        cum_V[6] * .kron_power(sigam, 6)
    ))
    cumulant_list[[6]] <- Cum6_GH
  }

  # Truncate list to r cumulants
  cumulant_list <- cumulant_list[1:r]

  if (nMu) cumulant_list <- Cum2Mom(cumulant_list, Type = "Multivariate")

  return(cumulant_list)
}


