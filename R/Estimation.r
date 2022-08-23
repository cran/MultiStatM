########### Estimation
## 0. Center_X - not visible in the help
## 1. conv_Stand_Multi
## 2. Esti_Skew_Mardia
## 3. Esti_Kurt_Mardia
## 4.  Esti_SkewVec_MRSz
## 5.  Esti_SkewInd_MRSz
## 6.  Esti_Skew_MRSz_2
## 7.  Esti_Kurt_CMRSz_2
## 8.  Esti_H3
## 9.  Esti_H4
## 10. Hermite_Nth
## 11. Skew_Mardia
## 12. Kurt_Mardia
## 13. Estimates_MMom_MCum
## 14 Esti_EVSK
## 15 Variance_of_Esti_Skew
## 16 Esti_Variance_Skew_Kurt
## 17 Esti_Hermite_Poly_HN_Multi
##


##  SkewEsti
##  SkewKron
##  SkewKronT
##  Variance_of
##################

# Centering a sample of vector variates,
#  @param x matrix of  sample, rows are observations of a d variate,
#  sample size is the number of rows
#  @return data matrix with rows centered by the sample means of columns
#  @family Standard
#
# @export
# Center_X <- function(x) {
#  apply(x, 2, function(y) y - mean(y))
  # this is: scale(x, center = TRUE, scale = FALSE)
#}


#' Standardize multivariate data
#'
#' For data formed by d-variate vectors x with sample covariance S and sample mean M,
#' it computes the values
#' \eqn{z=S^{-1/2}(x-M)}
#'
#' @param x a multivariate data matrix, sample size is the number of rows
#' @return a matrix of multivarate data with null mean vector and
#' identity sample covariance matrix
#' @examples
#' x<-MASS::mvrnorm(1000,c(0,0,1,3),diag(4))
#' z<-conv_Stand_Multi(x)
#' mu_z<- apply(z,2,mean)
#' cov_z<- cov(z)
#'
#' @export
conv_Stand_Multi<-function(x){
  # x is a multivariate vector (rows) of data, sample size is the number of rows,
  z<- scale(x, center = TRUE, scale = FALSE)  # Center_X(x)
  cx<-stats::cov(x)
  svdx<-svd(cx)
  sm12<-svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  z1<-t(sm12%*%t(z))
  return(z1)
}

##############################
#' Estimation of multivariate T-Moments and T-Cumulants
#'
#' Provides estimates of univariate and multivariate moments and cumulants up to order r.
#' By default data are standardized; using only demeaned or raw data is also possible.
#'
#' @param X d-vector data
#' @param r The highest moment order (r >2)
#' @param centering set to T (and scaling = F) if only centering is needed
#' @param scaling   set to T (and centering=F) if standardization of multivariate data is needed
#' 
#' @return \code{estMu.r}: the list of the multivariate moments up to order \code{r}  
#' @return \code{estCum.r}: the list of the multivariate cumulants up to order \code{r}
#' 
#' 
#' @family Estimation
#'
#' @examples
#' ## generate random data from a 3-variate skew normal distribution
#' alpha<-c(10,5,0)
#' omega<-diag(3)
#' x<-distr_SkewNorm_Rand(50,omega,alpha)
#' ## estimate the first three moments and cumulants from raw (uncentered and unstandardized) data
#' Esti_MMom_MCum(x,3,centering=FALSE,scaling=FALSE)
#' ## estimate the first three moments and cumulants from standardized data
#' Esti_MMom_MCum(x,3,centering=FALSE,scaling=TRUE)
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.
#' @export
Esti_MMom_MCum  <- function(X,r,centering = FALSE,scaling = TRUE){
  if (r<3) (stop("r must be greater than 2"))
  if (is.vector(X)) {
    m=mean(X)
    v=stats::var(X)
    if (centering==TRUE) X <- X-m  #Center_X(X)
    if (scaling == TRUE ) X <- (X-m)/sqrt(v)
    Mu_X_r <- NULL
    if (centering==T || scaling==T)  {Mu_X_r[1]<-0} else  {Mu_X_r[1] <- m}
    if (scaling==T) {Mu_X_r[2]<- 1} else  {Mu_X_r[2] <- mean(X^2)}
    for (j in c(3:r)) {
     Mu_X_r[j] <- mean(X^j)
    }
    Cum_X_r <- conv_Mom2Cum(Mu_X_r)
    EstiMomCum <- list(estMu.r=Mu_X_r,estCum.r=Cum_X_r)
    return(EstiMomCum)
     }
  if (dim(X)[1]<dim(X)[2]) (stop("X is not a proper data matrix"))
  d<-dim(X)[[2]]
  if (centering==TRUE) X <- scale(X, center = TRUE, scale = FALSE)  #Center_X(X)
  if (scaling == TRUE ) X <- conv_Stand_Multi(X)
  Mu_X_r <- NULL
  if (centering==TRUE || scaling==TRUE) {Mu_X_r[[1]]<-rep(0,d)} else  {Mu_X_r[[1]] <- colMeans(X)}
  for (j in c(2:r)) {
    Xadk <- apply(X, 1, function(y) .KronPower(y,j))
    if (scaling==T && j==2) { Mu_X_r[[j]]<-c(diag(d)) } else {Mu_X_r[[j]] <- rowMeans(Xadk)} #apply(Xadk, 2,mean)
  }
  Cum_X_r <- conv_Mom2CumMulti(Mu_X_r)
  EstiMomCum <- list(estMu.r=Mu_X_r,estCum.r=Cum_X_r)
  return(EstiMomCum)
}




####################
#' Estimation of multivariate Mean, Variance, T-Skewness and T-Kurtosis vectors
#'
#' Provides estimates of mean, variance, skewness and kurtosis vectors for  d-variate data
#' @param X d-variate data vector
#' @return The list of the estimated mean, variance, skewness and kurtosis vectors
#' @examples
#' x<- MASS::mvrnorm(100,rep(0,3), 3*diag(rep(1,3)))
#' EVSK<-Esti_EVSK(x)
#' names(EVSK)
#' EVSK$estSkew
#' @family Estimation
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, Sections 6.4.1 and 6.5.1
#' @export
Esti_EVSK <- function(X){
  if (is.vector(X)) {
    Mu_X<-mean(X)
    Vari_X<-stats::var(X)
    z<-(X-Mu_X)/sqrt(Vari_X)
    est.Skew<-mean(z^3)
    est.Kurt<-mean(z^4)
    estiEVSK <- list(Mu_X,Vari_X, est.Skew,est.Kurt )
    names(estiEVSK) <- c("estMu" ,   "estVar"  , "estSkew" , "estKurt")
    return(estiEVSK)
  }
  Mu_X <-  colMeans(X)  # column
  Vari_X <-stats::cov(X)
  kim <- Esti_MMom_MCum(X,4,centering=F,scaling = T)
  estiEVSK <- list(Mu_X,Vari_X, kim$estCum.r[[3]],kim$estCum.r[[4]] )
  names(estiEVSK) <- c("estMu" ,   "estVar"  , "estSkew" , "estKurt")
  return(estiEVSK)
}






#' Estimation of Mardia's Skewness index
#'
#' Compute the multivariate Mardia's skewness index and
#' provides the p-value for the hypothesis of zero symmetry under the
#' Gaussian assumption
#' @param x A matrix of multivariate data
#' @return \code{Mardia.Skewness} The skewness index
#' @return \code{p.value} The p-value  under the Gaussian hypothesis
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.1
#' @family Estimation
#' @export
Esti_Skew_Mardia<-function(x){
  z<-conv_Stand_Multi(x)
  n=dim(z)[1]
  d=dim(z)[2]
  MSkew<-sum((z%*%t(z))^3)/n^2
  pval<-stats::pchisq(n*MSkew/6,choose(d+2,3),lower.tail = FALSE)
  return(list("Mardia.Skewness"= MSkew,"p.value"= pval))
}


#' Estimation of Mardia's Kurtosis Index
#'
#' @param x A matrix of multivariate data
#' @return \code{Mardia.Kurtosis} The kurtosis index
#' @return \code{p.value} The p-value under the Gaussian hypothesis
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.1
Esti_Kurt_Mardia<-function(x){
  z<-conv_Stand_Multi(x)
  n <- dim(z)[1]
  d <- dim(x)[2]
  b2d <-sum(t(z^2)%*%z^2)/n
  std.norm <- sqrt(n)*(b2d-d*(d+2))/sqrt(8*d*(d+2))

  pval<- stats::pnorm(abs(std.norm), lower.tail = FALSE)*2
  return(list("Mardia.Kurtosis"= b2d,"p.value"= pval ))
}



#' Estimation of the Total Kurtosis Index
#'
#' @param x A matrix of multivariate data
#' @return \code{Total.Kurtosis} The total kurtosis index
#' @return \code{p.value} The p-value under the Gaussian hypothesis
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.1
Esti_Kurt_Total<-function(x){
  EVSK <- Esti_EVSK(x)
  MK <-sum(EVSK$estCurt^2)
  n <- dim(x)[1]
  d <- dim(x)[2]
  pval<-stats::pchisq(n*MK/24,choose(d+3,4),lower.tail = FALSE)
  return(list("Total.Kurtosis"= MK,"p.value"= pval ))
}



#' Estimation of Mori, Rohatgi, Szekely (MRSz's) skewness vector
#'
#' @param x A matrix of multivariate data
#' @return \code{MRSz.Skewness.Vector} The skewness vector
#' @return \code{MRSz.Skewness.Index} The skewness index
#' @return \code{p.value} The p-value for the hypothesis of symmetry under the
#' Gaussian assumption
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.2
#' @references  S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644
#' @export
Esti_Skew_MRSz<-function(x){
  z<-conv_Stand_Multi(x)
  n=dim(z)[1]
  d=dim(z)[2]
  z2<-apply(z^2,1,sum)
  MSv<-apply(z2*z,2,mean)
  MS<-sum(MSv^2)
  pval<-stats::pchisq(n*MS/(2*(d+2)),d,lower.tail = FALSE)
  return(list("MRSz.Skewness.Vector"=MSv,"MRSz.Skewness.Index"=MS,"p.value"=pval))
}

#' Estimation of Cardoso, Mori,  Rohatgi, Szekely (CMRSz's) kurtosis matrix
#'
#' @param x A matrix of multivariate data
#' @return \code{CMRSz.Kurtosis} The kurtosis matrix
#' @return \code{p.value} The p-value for the hypothesis of symmetry under the
#' Gaussian assumption
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.9
#'
#' @export
Esti_Kurt_CMRSz <- function(x){
   n=dim(x)[1]
  d<-dim(x)[2]

  Id <- diag(d)
  Id2 <- diag(d^2)
  szor <- kronecker(Id2,as.vector(Id))
  szort <- kronecker(Id2,t(as.vector(Id)))
  V4 <- szort %*%matr_Symmetry(d, 4)%*% szor
  V4.svd <- svd(V4)
  Pinv.V4.svdd<- c(1/V4.svd$d[1:(d*(d+1)/2) ], rep(0,d*(d-1)/2))
  P.inv <- V4.svd$u%*%diag(Pinv.V4.svdd)%*% t(V4.svd$v)/24
  evsk <- Esti_EVSK(x)
  M.Kurt <- as.vector (szort %*% evsk$estKurt)
  M.KurtInd <- t(M.Kurt)%*%P.inv%*%M.Kurt
  pval<-stats::pchisq(n*M.KurtInd,d*(d+1)/2,lower.tail = FALSE)
 return(list("CMRSz.Kurtosis"= matrix(M.Kurt,nrow = d),"p.value"=pval))
}



#' Asymptotic Variance of the estimated skewness vector
#'
#' @param cum The theoretical/estimated cumulants up to order 6 in vector form
#' @return The matrix of theoretical/estimated variance
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021. Ch.6, formula (6.13)
#' @examples
#' alpha<-c(10,5)
#' omega<-diag(rep(1,2))
#' MC <- distr_SkewNorm_MomCum_Th(r = 6,omega,alpha)
#' cum <- MC$CumX
#' VS <- Esti_Skew_Variance_Th(cum)
#' @family Estimation
#' @export
#'
Esti_Skew_Variance_Th <- function(cum){
  if (length(cum) < 6) (stop("cum  has  not large enough order"))
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j3
  d <- length(cum[[1]])
  type_2j3 <- c(0,0,2,0,0,0)
  # loc_type_2j3  <-  Partition_Type_eL_Location(type_2j3)
  L_2j3 <- .Commutator_Moment_eL(type_2j3,d)
  # %%%%%%%%%%%%%%%%%%%%%%%  Commutator  L_1j2_1j1
  # type_1j2_1j1  <- [1,1,0];%
  # loc_type_1j2_1j1   <-  PartitionType_eL_Location(type_1j2_1j1 );
  # L_1j2_1j1   <-  MomentCommutators_L(type_1j2_1j1 ,loc_type_1j2_1j1 (1),loc_type_1j2_1j1 (2),d);
  Id3 <-diag(d^3)
  K2 <-  kronecker( matr_Commutator_Kmn(d,d),diag(d)) #
  K3  <- matr_Commutator_Kmn(d,d^2);# L1_2_1_1 X van elõl
  Sym22 <-  Id3+K2+K3; #  Sym22 - L_1j2_1j1'
  #%%%%%%%%%%%%%%%% Commutator L2_H4

  L2_H4i  <-  kronecker(t(Sym22),t(Sym22))%*%matr_Commutator_Kperm(c(1,3,4,2,5,6 ),d)
  #L2_H4i  <-  kronecker(t(Sym22),t(Sym22))[indx_Commutator_Kperm(c(1,3,4,2,5,6 ),d)]

  #%%%%%%%%%%%%%%%5  Commutator M_3 inverse already!!!!
  d1 <- c(d,d,d)
  M3_m_ni  <- t(matr_Commutator_Mixing( d1,d1))
  #%%%%%%%%%%%%%%%%%%%%%
  Id <- diag(d)
  vec_Var_Skew  <- as.matrix(cum[[ 6]]) + t(L_2j3)%*%kronecker(as.matrix(cum[[ 3]]),as.matrix(cum[[ 3]])) +
    L2_H4i%*%kronecker(as.vector(Id),as.matrix(cum [[4]]))+
    M3_m_ni%*%kronecker(as.vector(Id), kronecker(as.vector(Id),as.vector(Id)))-
    kronecker(as.matrix(cum[[ 3]]),as.matrix(cum[[ 3]]))
  Var_Skew  <-  matrix(vec_Var_Skew, nrow=d^3)
  return(Var_Skew)
}







####################################

#' Estimated Variance of  skewness and kurtosis vectors
#'
#' Provides the estimated covariance matrices of the data-estimated skewness and kurtosis
#' vectors.
#'
#' @param X A matrix of d-variate data
#' @return The list of covariance matrices of the skewness and kurtosis vectors
# #' @examples
# #' d <- 2
# #' n <- 250
# #' x <-  matrix( 2*rnorm(d*n,3), nrow = n)
# #' evske <- Esti_Variance_Skew_Kurt(x)
# #' names(evske)
# #' vS <- matrix(evske$Vari_Skew_e,nrow=d^3)
# #' vK <- matrix(evske$Vari_Kurt_e,nrow=d^4)
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.
#' @export
Esti_Variance_Skew_Kurt<- function(X){
  estiEVSK <- Esti_EVSK(X) #
  Z <-  conv_Stand_Multi(X)
  H3Zt<- apply(Z,1,.Hermite_Third) # Hermite
  #
  cH3 <-  -  apply(H3Zt,2,function(U) U-estiEVSK$estSkew)
  # cov(t(H3Zt))
  Vari_Skew_e <- stats::cov(t(cH3))
  #############
  H4Zt <-  apply(Z,1,.Hermite_Fourth)
  # Est_H4(Z)
  cH4 <-  -  apply(H4Zt,2,function(U) U-estiEVSK$estCurt)
  Vari_Kurt_e <- stats::cov(t(cH4))
  esti.var.SK <- list(Vari_Skew_e, Vari_Kurt_e)
  names(esti.var.SK) <- c("Vari_Skew_e" ,   "Vari_Kurt_e")
  return(esti.var.SK)
}

###################

#' Asymptotic Variance of the estimated  kurtosis vector
#'
#' @param cum The theoretical/estimated cumulants up to the 8th order in vector form
#' @return The matrix of theoretical/estimated variance
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. Ch. 6, formula (6.26)
#' @family Estimation
#' @export

Esti_Kurt_Variance_Th <- function(cum){
  if (length(cum) < 8) (stop("cum  has  not large enough order"))
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_1j3_1j5
  d <- length(cum[[1]])
  # kron2 <- function(Mm) kronecker(Mm,Mm)
  Id <- diag(d)
  vId <- as.vector(Id)
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j3
  type_2j3 <- c(0,0,2,0,0,0)
 # loc_type_2j3  <-  Partition_Type_eL_Location(type_2j3)
  L_2j3 <- .Commutator_Moment_eL(type_2j3,d)
  ###########L_1j3_1j5
  type_1j3_1j5 <- c(0,0,1,0,1,0,0,0)
 # loc_type_1j3_1j5  <-  Partition_Type_eL_Location(type_1j3_1j5)
  L_1j3_1j5 <- .Commutator_Moment_eL(type_1j3_1j5,d)
  # %%%%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j4
  type_2j4 <- c(0,0,0,2,0,0,0,0)
  # loc_type_2j4  <-  Partition_Type_eL_Location(type_2j4)
  L_2j4 <- .Commutator_Moment_eL(type_2j4,d)
  #######
  ############## L2_H6i

  L_1j1_1j3 <-  diag(d^4) + matr_Commutator_Kperm(c(2, 1, 3, 4),d)+
    matr_Commutator_Kperm(c(3, 1, 2, 4) ,d)+
    matr_Commutator_Kperm(c( 4, 1, 2, 3),d)
  Legy <- kronecker(L_1j1_1j3,L_1j1_1j3)

 L2_H6i=Legy*matr_Commutator_Kperm(c(1, 3, 4,  5, 2, 6, 7, 8),d)
   # L2_H6i=Legy[indx_Commutator_Kperm(c(1, 3, 4,  5, 2, 6, 7, 8),d)]

  #%%%%%%%%%%%%%%%5  Commutator M_4_m_n
  d1 <- rep(d,4)
  M4_m_ni  <- t(matr_Commutator_Mixing( d1,d1))
  #%%%%%%%%%%%%%%%%%%%%% L_1j3_1j1
  type_1j3_1j1=c(1,0,1,0)
  # loc_type_1j3_1j1  <- Partition_Type_eL_Location(type_1j3_1j1)
  L_1j3_1j1=.Commutator_Moment_eL(type_1j3_1j1,d)
  ################
  L22_H4m<- .L22_H4(d)

  krId2 <- .kron2(vId)
  # %% kurtosis (6.26)
  vec_Var_Kurt =cum[[ 8]] +
    t(L_1j3_1j5)%*%kronecker(cum[[5]],cum[[3]])+
    t(L_2j4)%*%.kron2(cum[[4]]) - as.vector(.kron2(cum[[4]])) +
    L2_H6i%*%kronecker(vId,cum[[ 6]]+
                         t(L_2j3)%*%.kron2(cum[[3]]))+
    t(L22_H4m)%*%kronecker(krId2,cum[[4]])+
    M4_m_ni%*%.kron2(krId2)+
    .kron2(t(L_1j3_1j1))%*%kronecker(.kron2(cum[[3]]),vId)

  Var_Kurt=matrix(vec_Var_Kurt, nrow= d^4);
  return(Var_Kurt)
}


#' Estimate the N-th d-variate Hermite polynomial
#'
#' The vector x is standardized and the N-th d-variate polynomial is computed
#'
#' @param x a d-variate data vector
#' @param N the order of the d-variate Hermite polynomial
#' @return The vector of the N-th d-variate polynomial
#' @examples
#' x<-MASS::mvrnorm(100,rep(0,3),diag(3))
#' H3<-Esti_Hermite_Poly_HN_Multi(x,3)
Esti_Hermite_Poly_HN_Multi<-function(x,N){
  z<-conv_Stand_Multi(x)
  HN<-apply(apply(z,1,Hermite_Nth, N=N),1,mean)
  return(HN)
}
