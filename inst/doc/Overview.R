## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MultiStatM)

## -----------------------------------------------------------------------------
PTA<-Partition_Type_All(4)

## -----------------------------------------------------------------------------
PTA$Part.class[[1]]
PTA$Part.class[[2]]
PTA$Part.class[[3]] 

## -----------------------------------------------------------------------------
PTA$S_N_r 
PTA$eL_r

## -----------------------------------------------------------------------------
PTA$S_r_j[[2]] ## Partitions with size r=2, includes two types (above) each with number

## -----------------------------------------------------------------------------
a1<-c(1,2)
a2<-c(2,3,4)
a3<-c(1,3)
p1<-a1%x%a2%x%a3
as.vector(matr_Commutator_Kperm(c(3,1,2),c(2,3,2))%*%p1) ## same result as below
a3%x%a1%x%a2   

## -----------------------------------------------------------------------------
p1[indx_Commutator_Kperm(c(3,1,2),c(2,3,2))]

## -----------------------------------------------------------------------------
a<-c(1,2)
a3<-a%x%a%x%a
a3
as.vector(matr_Elimination(2,3)%*%a3)
as.vector(matr_Qplication(2,3)%*%matr_Elimination(2,3)%*%a3)

## -----------------------------------------------------------------------------
x<-c(1,2,3)
H2<-Hermite_Poly_HN_Multi(x,N=2)
H2[[1]]
H2[[2]]

## -----------------------------------------------------------------------------
H2<-Hermite_Poly_HN_Multi(x,Sig2=4*diag(3),N=2)
H2[[1]]
H2[[2]]

## -----------------------------------------------------------------------------

Hermite_Poly_NH_Multi_Inv(H2,N=2,Sig2=4*diag(3))[[1]] 

## -----------------------------------------------------------------------------
Covmat<-matrix(c(1,0.8,0.3,0.8,2,1,0.3,1,2),3,3)
Cov_X1_X2 <- Hermite_N_Cov_X1_X2(SigX12=Covmat,N=3)

## -----------------------------------------------------------------------------
mu<-list(c(1,1),c(2,1.5,1.5,2),c(4,3,3,3,3,3,3,4),c(10,7,7,6.5,7,6.5,6.5,7,7,6.5,6.5,7,6.5,7,7,10))
cum<-conv_Mom2CumMulti(mu)
cum

## -----------------------------------------------------------------------------
conv_Cum2MomMulti(cum)

## -----------------------------------------------------------------------------
mu[[3]][indx_Elimination(2,3)]

## -----------------------------------------------------------------------------
r.mu<-matr_Elimination(2,3)%*% mu[[3]]
as.vector(r.mu)

## -----------------------------------------------------------------------------
as.vector(matr_Qplication(2,3)%*%r.mu)

## -----------------------------------------------------------------------------
r.mu[indx_Qplication(2,3)]

## -----------------------------------------------------------------------------
alpha<-c(10,5,0)
omega<-diag(3)
MSN<-distr_SkewNorm_MomCum_Th(r=3,omega,alpha,nMu=TRUE)
round(MSN$Mu[[3]],3)
round(MSN$CumX[[3]],3)

## -----------------------------------------------------------------------------
distr_Uni_EVSK_Th(3, nCum = TRUE)$Kurt.U

## -----------------------------------------------------------------------------
data<-distr_SkewNorm_Rand(1000,omega,alpha)
EsMSN<-Esti_EVSK(data)
ThMSN<-distr_SkewNorm_EVSK_Th(omega,alpha)

## -----------------------------------------------------------------------------
EsMSN$estSkew[indx_Elimination(3,3)]   
ThMSN$SkewX[indx_Elimination(3,3)]   


## -----------------------------------------------------------------------------
EsMSN$estSkew[indx_UnivMomCum(3,3)]  ## Get univariate skewness for X1,X2,X3
EsMSN$estKurt[indx_UnivMomCum(3,4)]  ## Get univariate kurtosis for X1,X2,X3

## -----------------------------------------------------------------------------
Esti_Skew_Mardia(data)
ThMSN$SkewX.tot

## -----------------------------------------------------------------------------
Esti_Skew_MRSz(data)

## -----------------------------------------------------------------------------
as.vector(t(c(diag(3))%x%diag(3))%*%ThMSN$SkewX)  ## Theoretical MRS skewness vector

