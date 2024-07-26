########### Hermite
## 1. Hermite_Coeff
## 2. Hermite_Poly_HN
## 3. Hermite_CoeffMulti
## 4. Hermite_Poly_HN_Multi
## 5. Hermite_Poly_NH_Inv
## 6. Hermite_Poly_NH_Multi_Inv
## 7. Hermite_n_Cov_X1_X2






.Hermite_Coeff<- function(N){
  Hcoeff<-rep(0,floor(N/2)) #ceiling(N/2+1-(N%%2))
  for (k in c(0:floor(N/2))) {
    Hcoeff[k+1] <- (-1)^k*factorial(N)/factorial(N-2*k)/factorial(k)/2^k
  }

  return(Hcoeff)
}



.Hermite_Poly_HN<-function(x,N,sigma2=1){
  H_N_x<-rep(0,N)
  H_N_x[1]<-x
  for (n in 2:N){
    Hcoeff<-.Hermite_Coeff(n)
    nH<-length(Hcoeff)
    powersX<-seq(n,0,by=-2)
    powersSigma2<-(n-powersX)/2
    x_powers=rep(x,nH)^powersX
    sigma2_powers<-rep(sigma2,nH)^powersSigma2
    H_N_x[n]<-sum(Hcoeff*x_powers*sigma2_powers)
  }
  return(H_N_x)
}



.Hermite_Poly_NH_Inv<-function(H_N_x,sigma2=1){
  N<-length(H_N_x)
  x_val<-rep(0,N)
  x_val[1]=H_N_x[1]
  for (n in 2:N) {
    Hcoeff<-.Hermite_Coeff(n)
    nH<-length(Hcoeff)
    powersX<-seq(n,0,by=-2)
    if ((n%%2==0)) {H_N_x1 =rev(c(1,H_N_x[seq(2,n,by=2)]))}
    else {H_N_x1 =rev(H_N_x[seq(1,n,by=2)])}

    powersSigma2<-(n-powersX)/2
    signCoeff<-rep(-1,nH)^powersSigma2
    Xcoeff<-Hcoeff*signCoeff
    sigma2_powers<-rep(sigma2,nH)^powersSigma2
    x_val[n]<-sum(Xcoeff*H_N_x1*sigma2_powers)
  }
  return(x_val)
}


.Hermite_CoeffMulti<-function(N,d){
  PTA<-PartitionTypeAll(N)
  el_j<-PTA$eL_r
  HcoeffMatrix<-vector(mode = "list", length = ceiling(N/2+1-(N%%2)))
  kk=0
  for (k in 0:N) {
    if (N%%2== k%%2) {
      el=c(k,(N-k)/2,rep(0,N-2))
      loc_type_el<-c(0,0)
      for (m in 1:N) {
        if (is.vector(el_j[[m]])){
          if (prod((el==el_j[[m]]))) {loc_type_el<-c(m,1)}
        }
        else {
          for (mm in 1:dim(el_j[[m]])[1])
            if (prod((el==el_j[[m]][mm,]))) {loc_type_el<-c(m,mm)}
        }

      }
      kk=kk+1;
      HcoeffMatrix[[kk]]<- (-1)^((N-k)/2)*t(.Commutator_Moment_eL(el,d))
    }
  }
  HcoeffMatrix<-HcoeffMatrix[seq(ceiling(N/2+1-(N%%2)),1,by=-1)]
  return(HcoeffMatrix)
}


.Hermite_Poly_HN_Multi<-function(x,N,Sig2=diag(length(x))){

  d=length(x)
  H_N_x<-vector(mode = "list", length = N)
  x_ad<- vector(mode="list",length=N)
  x_ad[[1]]=x
  if (N>1){
    for (n in 2:N){
      x_ad[[n]]<-kronecker(x_ad[[n-1]],x)
    }
  }

  vSig2<-c(Sig2)
  vSig2_ad<-vector(mode="list",length=ceiling(N/2+1-(N%%2)))
  vSig2_ad[[1]]<-1
  if (N>1){
    for (n in 2:ceiling(N/2+1-(N%%2))){
      vSig2_ad[[n]]<-kronecker(vSig2_ad[[n-1]],vSig2)
    }
  }

  H_N_x[1]<-x_ad[1]
  if (N>1) {
    for (n in 2:N){
      HcoeffMatrix<-.Hermite_CoeffMulti(n,d)
      nH=length(HcoeffMatrix)
      X_powers<-vector(mode="list",length=nH)

      if ((n%%2)==0){
        X_powers[1:(nH-1)]<-x_ad[seq(n,1,by=-2)]
        X_powers[[nH]]<-1
      }
      else {X_powers<-x_ad[seq(n,1,by=-2)]}

      Sigma2_powers<-vSig2_ad[1:nH]
      H_N_x0=0
      for (k in 1:nH){
        H_N_x0<-H_N_x0+HcoeffMatrix[[k]]%*%kronecker(Sigma2_powers[[k]],X_powers[[k]])
      }
      H_N_x[[n]]=as.vector(H_N_x0)
    }
  }

  return(H_N_x)
}


.Hermite_Poly_NH_Multi_Inv<-function(H_N_X,N,Sig2=diag(length(H_N_X[[1]]))) {

  d<-length(H_N_X[[1]])
  x_ad_val<-vector(mode="list",length=N)
  vSig2<-c(Sig2)
  vSig2_ad<-vector(mode="list",length = ceiling(N/2+1-(N%%2)))
  vSig2_ad[[1]]<-1
  if (N>1) {
    for (n in 2:ceiling(N/2+1-(N%%2))){
      vSig2_ad[[n]]<-kronecker(vSig2_ad[[n-1]],vSig2)
    }
  }
  x_ad_val[1]<-H_N_X[1]

  if (N>1){
    for (n in 2:N) {
      HcoeffMatrix<-.Hermite_CoeffMulti(n,d)
      nH<-length(HcoeffMatrix)
      H_N_Xp<-vector(mode="list",length=nH)
      if ((n%%2==0)){
        H_N_Xp[1:(nH-1)]<-H_N_X[seq(n,1,by=-2)]
        H_N_Xp[[nH]]<-1
      }
      else {H_N_Xp <- H_N_X[seq(n,1,by=-2)]}
      Sigma2_powers<-vSig2_ad[1:nH]
      X0<-0
      for (k in 1:nH) {
        X0<-X0+(-1)^(k-1)*HcoeffMatrix[[k]]%*%kronecker(Sigma2_powers[[k]],H_N_Xp[[k]])
      }
      x_ad_val[[n]]<-as.vector(X0)
    }
  }
  return(x_ad_val)

}


#' Covariance matrix  for multivariate  T-Hermite polynomials
#'
#' Computation of the covariance matrix between d-variate T-Hermite polynomials
#' \eqn{H_N(X_1)} and \eqn{H_N(X_2)}.
#' @param SigX12 Covariance matrix  of the Gaussian vectors X1 and X2 respectively
#' of dimensions d1 and  d2
#' @param N Common degree of the multivariate Hermite polynomials
#' @return Covariance matrix of \eqn{H_N(X_1)} and \eqn{H_N(X_2)}
#'
#' @examples
#' Covmat<-matrix(c(1,0.8,0.8,1),2,2)
#' Cov_X1_X2 <- HermiteCov12(Covmat,3)
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. (4.59),  (4.66),
#'
#' @family Hermite Polynomials
#' @export

HermiteCov12 <- function(SigX12,N){
  #
  dimX <- dim(SigX12)
  d1 <- rep(dimX[1],N)
  d2 <- rep(dimX[2],N)

  vSig2 <- as.vector(SigX12)
  vSig22 <- .KronPower(vSig2,N);
  vSig23 <-  .indx_Commutator_Mixing_t(vSig22,d1,d2)
  CH_1_2_n  <-  matrix(vSig23,nrow=d1[1]^N )
  return(CH_1_2_n)
}


.Hermite_Nth <-function(x,N){
  d=length(x)
  HN<-.Hermite_Poly_HN_Multi(x,N,diag(d))[[N]]
  return(as.vector(HN))
}





