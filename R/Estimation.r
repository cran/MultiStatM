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
#  @param x matrix of  sample, rows are observations of a d variate random vector,
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
#' @return a matrix of multivariate data with null mean vector and
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

#' Estimation of  Mori,  Rohatgi, Szekely (CMRSz's) kurtosis matrix
#'
#' @param x A matrix of multivariate data
#' @return \code{MRSz.Kurtosis} The kurtosis matrix
#' @return \code{p.value} The p-value for the hypothesis of symmetry under the
#' Gaussian assumption
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.9
#'
#' @export
Esti_Kurt_MRSz <- function(x){
  n=dim(x)[1]
  d<-dim(x)[2]
  Id <- diag(d)
  Id2 <- diag(d^2)
  i_ind <- NULL 
  j_ind <-  NULL 
  col_ind <-  NULL
  for (k in 0: (d-2)){
    col_ind <- c(col_ind, (ceiling(c((k*(d+1)+2):((k+1)*(d+1)))/d)-1)*d + 
                   ((k*(d+1)+2):((k+1)*(d+1)) %% d)  + 
                   ((c((k*(d+1)+2):((k+1)*(d+1))) %% d)==0)*d)
    i_ind <- c(i_ind, (ceiling(c((k*(d+1)+2):((k+1)*(d+1)))/d)))
    j_ind <- c(j_ind, ((k*(d+1)+2):((k+1)*(d+1)) %% d)  + 
                 ((c((k*(d+1)+2):((k+1)*(d+1))) %% d)==0)*d)
  }
  Matr <- matrix(rep(0,d^4),nrow = d^2)
  Sz_all <- c(Id)
  for (i in c(0:(d-1))){
    Matr[ ,i*(d+1)+1] <- (4*(d-1)+24)*kronecker(Id[ ,i+1],Id[ ,i+1])+
      4*(Sz_all-kronecker(Id[ ,i+1],Id[ ,i+1]));
  }
  sz <- 1
  for (kk in col_ind){
    
    Matr[ ,kk] <- (2*(d-2)+12)*(kronecker(Id[ ,i_ind[sz]],Id[ ,j_ind[sz]])+
                                  kronecker(Id[,j_ind[sz]],Id[,i_ind[sz]]))
    sz <- sz+1
  }
  
  szor <- kronecker(Id2,as.vector(Id))
  szort <- kronecker(Id2,t(as.vector(Id)))
  V4 <- Matr #/24
  V4.svd <- svd(V4)
  si <- sqrt( 1/V4.svd$d[1:(d*(d+1)/2) ])
  Vi <- V4.svd$v[,1:(d*(d+1)/2)]
  #  Pinv.V4.svdd<- c(, rep(0,d*(d-1)/2))
  #  P.inv <- V4.svd$u%*%diag(Pinv.V4.svdd)%*% t(V4.svd$v)/24
  evsk <- Esti_EVSK(x)
  M.Kurt <- as.vector (szort %*% evsk$estKurt)
  M.Kurti <- diag(si)%*%t(Vi)%*%M.Kurt
  M.KurtInd <- t(M.Kurti)%*%M.Kurti  #/24
  #t(M.Kurt)%*%P.inv%*%M.Kurt
  pval<-stats::pchisq(n*M.KurtInd,d*(d+1)/2,lower.tail = FALSE)
  return(list("MRSz.Kurtosis"= matrix(M.Kurt,nrow = d),"pval"=pval))
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
  d <- length(cum[[1]])
  type_2j3 <- c(0,0,2,0,0,0)
  Id <- diag(d)
  perm1 <- c(2,1,3)
  perm2 <- c(2,3,1) # inverse
  D <- d^(2*length(perm1))
  
  ########  FIRST TERM ###########
  B1<-as.vector(.indx_Commutator_Moment_t(kronecker(cum[[ 3]],cum[[ 3]]),type_2j3,d)) 
  
  #### SECOND TERM###############
  Y <- kronecker(as.vector(diag(d)),as.vector(cum [[4]]))[indx_Commutator_Kperm(c(1,3,4,2,5,6 ),d)]
  B2 <- rep(0,D)
  permM <- matrix(c(1:3,perm1,perm2),nrow=3,byrow=TRUE)
  for(k in 1:3){
    for (m in 1:3) {
      B2 <- B2 + Y[.indx_kron2(permM[k,],permM[m,],d)] 
    }
  }
  
  ########################### THIRD TERM #############
  d1 <- c(d,d,d)
  xp<-kronecker(as.vector(Id), kronecker(as.vector(Id),as.vector(Id)))
  B3<-xp[indx_Commutator_Kperm(c(1,3,5,2,4,6 ),d)]+xp[indx_Commutator_Kperm(c(1,3,5,2,6,4 ),d)]+xp[indx_Commutator_Kperm(c(1,3,5,4,2,6 ),d)]+
    xp[indx_Commutator_Kperm(c(1,3,5,6,2,4 ),d)]+xp[indx_Commutator_Kperm(c(1,3,5,4,6,2 ),d)]+xp[indx_Commutator_Kperm(c(1,3,5,6,4,2 ),d)]
  
  #### Variance formula 6.13
  vec_Var_Skew  <- cum[[6]] + B1 + B2 + B3 - kronecker(cum[[ 3]],cum[[ 3]])
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
#' Warning: the function requires 8! computations, for d>3, the timing required maybe large.  
#' 
#' @param cum The theoretical/estimated cumulants up to the 8th order in vector form
#' @return The matrix of theoretical/estimated variance
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. Ch. 6, formula (6.26)
#' @family Estimation
#' @export

Esti_Kurt_Variance_Th <- function(cum){
  if (length(cum) < 8) (stop("cum  has  not large enough order"))
  d <- length(cum[[1]])
  D <- d^(2*4)
  Id <- diag(d)
  vId <- as.vector(Id)
  krId2 <- kronecker(vId,vId)
  
  ########################### First term
  B1 <- cum[[ 8]] 
  
  #%%%%%%%%%%%%%%%%%%%% Second  Commutator  L_1j3_1j5
  type_1j3_1j5 <- c(0,0,1,0,1,0,0,0)
  B2 <- indx_Commutator_Moment(kronecker(cum[[5]],cum[[3]]),
                               type_1j3_1j5,d)
  
  # %%%%%%%%%%%%%%%%%%%%%%% Third  Commutator  L_2j4 
  type_2j4 <- c(0,0,0,2,0,0,0,0)
  B3 <- .indx_Commutator_Moment_t(kronecker(cum[[4]],cum[[4]]),
                                  type_2j4,d)
  
  #%%%%%%%####################### Fourth  Commutator  L_2j3, L2_H6i
  type_2j3 <- c(0,0,2,0,0,0)
  x0 <- cum[[6]] + .indx_Commutator_Moment_t(kronecker(cum[[3]],cum[[3]]),
                                             type_2j3,d)
  x1 <- kronecker(vId,x0)
  x <- x1[indx_Commutator_Kperm(c(1, 3, 4,  5, 2, 6, 7, 8),d)]
  # see (A.13) for L_2H6
  permM <- matrix(c(1:4,c(2, 1, 3, 4),c(3, 1, 2, 4),c( 4, 1, 2, 3)),
                  nrow=4,byrow=TRUE)
  B4 <- rep(0,D)
  for(k in 1:4){
    for (m in 1:4) {
      B4 <- B4 + x[.indx_kron2(permM[k,],permM[m,],d)] 
    }
  }
  ############## Fifth
  x <-  kronecker(krId2,cum[[4]])
  B5 <-  .indx.L22_H4_t(x,d)
  
  #%%%%%%%%%%%%%%%########### Sixth   Commutator M_4_m_n
  d1 <- rep(d,4)
  x <- kronecker(krId2,krId2)
  B6  <-.indx_Commutator_Mixing_t(x,d1,d1)
  
  #%%%%%%%%%%%%%####### Seven 
  B7  <- - kronecker(cum[[4]],cum[[4]])
  
  ###################################### Eith  %% L_1j3_1j1
  type_1j3_1j1=c(1,0,1,0)
  x <- kronecker(kronecker(cum[[3]],cum[[3]]),vId)
  permM <- matrix(c(1:4,c(1,2,4,3),c(1,4,2,3),c(4,1,2,3)),
                  nrow=4,byrow=TRUE) # inverse
  B8 <- rep(0,D)
  for(k in 1:4){
    for (m in 1:4) {
      B8 <- B8 + x[.indx_kron2(permM[k,],permM[m,],d)]
    }
  }
  ################
  # %% kurtosis (6.26)
  vec_Var_Kurt  <-  B1+B2+B3+B4+B5+B6+B7+B8 
  Var_Kurt <- matrix(vec_Var_Kurt, nrow= d^4);
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

#' Gram-Charlier approximation to a multivariate density
#' 
#' Provides the truncated Gram-Charlier approximation to a multivariate density. Approximation can
#' be up to the first k=8 cumulants. 
#' 
#' @param X A matrix of d-variate data
#' @param k the order of the approximation, by default set to 4; 
#' (k must not be smaller than 3 or greater than 8)
#' @param cum if NULL (default) the cumulant vector is estimated from X. 
#' If \code{cum} is provided no estimation of cumulants is performed.
#' @return The vector of the Gram-Charlier density evaluated at X
#' 
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021. Section 4.7.
#' 
#' @examples
#' # Gram-Charlier density approximation (k=4) of data generated from 
#' # a bivariate skew-gaussian distribution
#' n<-50
#' alpha<-c(10,0)
#' omega<-diag(2)
#' X<-distr_SkewNorm_Rand(n,omega,alpha)
#' EC<-Esti_EVSK(X)
#' fy4<-Esti_Gram_Charlier(X[1:5,],cum=EC)
#' @export
Esti_Gram_Charlier<-function(X,k=4,cum=NULL){
  if (!is.null(cum)) {k=length(cum)}
  if (k<3) stop("k must be greater than 2")
  if (k>8) stop("k cannot be greater than 8")
  if (is.vector(X)) stop(" X must be a data matrix")
  
  d<-dim(X)[[2]]
  
  if (!is.null(cum)) {EC<-cum
  z1<-t(apply(X,1, function(x) x-as.vector(EC[[1]]))) 
  if (is.vector(EC[[2]])) {cx<-matrix(EC[[2]],nrow=d)} else {cx<-EC[[2]]}
  svdx<-svd(cx)
  sm12<-svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  if (is.vector(z1)) {as.matrix(z1);Z<-t(sm12%*%z1) }  else {Z<-t(sm12%*%t(z1))}
  } 
  else {EC<-Esti_MMom_MCum(X,k)$estCum.r
  Z<-conv_Stand_Multi(X)
  } 
  
  gy<-1
  for (j in 3:min(k,5)){
    HN<-apply(Z,1,Hermite_Nth, N=j)
    gy<-gy+EC[[j]]%*%HN/factorial(j)
  }
  phi<-mvtnorm::dmvnorm(Z,rep(0,d),diag(d))
  
  if (k>5){
    if (k==6) {B6<-indx_Symmetry(EC[[6]]+10*kronecker(EC[[3]],EC[[3]]),d,6)
    Bell<-list(B6)
    }
    if (k==7) {B6<-indx_Symmetry(EC[[6]]+10*kronecker(EC[[3]],EC[[3]]),d,6)
    B7<-indx_Symmetry(EC[[7]]+35*kronecker(EC[[3]],EC[[4]]),d,7)
    Bell<-list(B6,B7)
    }
    if (k==8) {B6<-indx_Symmetry(EC[[6]]+10*kronecker(EC[[3]],EC[[3]]),d,6)
    B7<-indx_Symmetry(EC[[7]]+35*kronecker(EC[[3]],EC[[4]]),d,7)
    B8<-indx_Symmetry(EC[[8]]+56*kronecker(EC[[3]],EC[[5]])+35*kronecker(EC[[4]],EC[[4]]),d,8)
    Bell<-list(B6,B7,B8)
    }
    for (j in 6:min(k,8)){
      HN<-apply(Z,1,Hermite_Nth, N=j)
      gy<-gy+Bell[[j-5]]%*%HN/factorial(j)
    }
  }
  
  GC<-as.vector(gy*phi)
  return(GC)
}

