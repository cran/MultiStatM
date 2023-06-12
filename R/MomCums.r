##### Moments and cumulants
## 1. conv_Mom2Cum
## 2. conv_Cum2Mom
## 3. conv_Mom2CumMulti
## 4. conv_Cum2MomMulti
##


#' Convert moments to cumulants (univariate)
#'
#' Obtains a vector of univariate cumulants from a vector of
#' univariate moments
#'
#' @param mu_x the r-vector of moments starting from the first moment - the mean -
#' and arriving to the r-th order moment
#'
#' @return \code{cum_x} the vector of univariate cumulants
#'
#' @examples
#' mu_x<- c(1,2,3,4)
#' conv_Mom2Cum(mu_x)
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section 3.4 formula 3.20
#'
#' @family Moments and cumulants
#'
#' @export
conv_Mom2Cum<-function(mu_x){
N<-length(mu_x)
cum_x<-rep(0,N)
for (n in 1:N){
  PT<-Partition_Type_All(n)
  eL_r<-PT$eL_r
  S_r_j<-PT$S_r_j
  cum=0
  for (m in 1:n) {
    if (is.vector(eL_r[[m]])) {
      mu=mu_x[1:length(eL_r[[m]])]
      cum<- cum +(-1)^(m-1)*factorial(m-1)*S_r_j[[m]][1]*prod((mu[eL_r[[m]]!=0])^eL_r[[m]][(eL_r[[m]]!=0)])}
   else  {
     mk<-min(dim(eL_r[[m]]))
     mu<-mu_x[1:ncol(eL_r[[m]])]
     for (k in 1:mk){
           cum<-cum+(-1)^(m-1)*factorial(m-1)*S_r_j[[m]][k]*prod((mu[eL_r[[m]][k,]!=0])^eL_r[[m]][k,(eL_r[[m]][k,]!=0)])
     }
   }
  }
  cum_x[n]<-cum
}
return("Cum"=cum_x)
}


#' Convert cumulants to moments (univariate)
#'
#' Obtains a vector of univariate moments from a vector of
#' univariate cumulants
#'
#' @param cum_x the r-vector of cumulants starting from the first - the mean -
#' and arriving to the r-th order cumulant
#'
#' @return \code{mu_x} the vector of univariate moments
#'
#' @examples
#' cum_x<- c(1,2,3,4)
#' conv_Cum2Mom(cum_x)
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section 3.4 formula 3.23
#'
#' @family Moments and cumulants
#'
#' @export
conv_Cum2Mom <-function(cum_x){
  N<-length(cum_x)
  mu_x<-rep(0,N)
  for (n in 1:N){
    PT<-Partition_Type_All(n)
    eL_r<-PT$eL_r
    S_r_j<-PT$S_r_j
    mu=0
    for (m in 1:n) {
      if (is.vector(eL_r[[m]])) {
        cum=cum_x[1:length(eL_r[[m]])]
        mu<- mu +S_r_j[[m]][1]*prod((cum[eL_r[[m]]!=0])^eL_r[[m]][(eL_r[[m]]!=0)])}
      else  {
        mk<-min(dim(eL_r[[m]]))
        cum<-cum_x[1:ncol(eL_r[[m]])]
        for (k in 1:mk){
          mu<-mu+S_r_j[[m]][k]*prod((cum[eL_r[[m]][k,]!=0])^eL_r[[m]][k,(eL_r[[m]][k,]!=0)])
        }
      }
    }
    mu_x[n]<-mu
  }
  return("Mom"=mu_x)
}





#' Convert T-moments to T-cumulants  (multivariate)
#'
#' Obtains a vector of d-variate cumulants from a vector of
#' d-variate moments
#'
#'
#' @param mu the list of r d-variate moments in vector form starting from the first moment -
#' the vector of means - and arriving to the r-th order d-variate moment in vector form
#'
#' @return \code{Cum} the list of n vectors of d-variate cumulants
#'
#' @examples
#' #Mu contains the  T-vector moments up to the fifth order of the bivariate
#' #standard normal distribution
#' mu<-list(c(0,0),c(1,0,0,1),c(0,0,0,0,0,0,0,0),c(3,0,0,1,0,1,1,0,0,1,1,0,1,0,0,3),
#' c(rep(0,32)))
#' conv_Mom2CumMulti(mu)
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section 3.4 formula 3.29
#'
#' @family Moments and cumulants
#'
#' @export
#'
conv_Mom2CumMulti = function(mu){
  N  <-  length(mu)
  d <- length(mu[[1]])

  cum  <-  vector(mode = "list", length = N)
  cum[[1]] <-  mu[[1]]

  for(n in 2:N){
    PA <- Partition_Type_All(n)
    el_js <-  PA$eL_r
    # eL_r<-PT$eL_r
    S_r_j<-PA$S_r_j
    an <-   1:n
    L_szor_mu  <-  array(0,c(d^n,n))
    L_szor_mu[,1] <- mu[[n]]

    if ( n>=3 ) {
      for(r in 2:(n-1)){

        if(is.null(dim(el_js[[r]]))){
          el_jsr <- t(as.matrix(el_js[[r]], c(1,length(el_js[[r]]))))### vettore riga ora
          #S_r_jr <-as.vector(S_r_j[[r]])
        }else{
          el_jsr <- el_js[[r]]
          # S_r_jr <- S_r_j[[r]]
        }
        Prod_M_k <- array(0, c(d^n, dim(el_jsr)[1]))
        for(m in 1:dim(el_jsr)[1]){
          el <- el_jsr[m,]
          j_el <- an[el!=0]
          el_j <- el[el!=0]
          Kr_1 <- vector(mode = "list", length = length(j_el))
          for(k in 1:length(j_el)){
            Mc1 <-array(mu[j_el[k]], c(el_j[k],el_j[k])) 
            Kr_1[[length(j_el)-k+1]] <- .KronProd(Mc1[1,])
          }
          Kr_m <- S_r_j[[r]][m]* .KronProd(Kr_1)
          Prod_M_k[,m] <- as.matrix(Kr_m)# L_el%*%
        }

        L_szor_mu[,r]<-(-1)^(r-1)*factorial(r-1)*apply(Prod_M_k,1,sum)
      }
    }
    McN <- mu[rep(1,n)]
    L_szor_mu[,n] <- (-1)^(n-1)*factorial(n-1)*as.matrix(.KronProd(McN))
    # S_r_j[[n]][1]=1  !!!
    #cum[[n]] <- as.vector(t(as.matrix(MSymm%*% as.matrix(apply(L_szor_mu,1,sum)))))
    xyz<-apply(L_szor_mu,1,sum)
    cum[[n]]<-indx_Symmetry(xyz,d,n)
  }

  return("Cum"= cum)
}






#' Convert T-cumulants to T-moments   (multivariate)
#'
#' Obtains a vector of d-variate moments from a vector of
#' d-variate cumulants
#'
#'
#' @param cum the list of r d-variate cumulants in vector form starting from the first cumulant -
#' the vector of means - and arriving to the r-th order d-variate cumulant in vector form
#'
#' @return \code{Mom}  the list of r vectors of d-variate moments
#'
#' @examples
#' #cum contains the  T-vector cumulants up to the fifth order of the bivariate
#' #standard normal distribution
#' cum<-list(c(0,0),c(1,0,0,1),c(0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#' c(rep(0,32)))
#' conv_Cum2MomMulti(cum)
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Section 3.4 formula 3.40
#'
#' @family Moments and cumulants
#'
#' @export
#'
conv_Cum2MomMulti = function(cum){
  N = length(cum)
  d = length(cum[[1]])

  mu = vector(mode = "list", length = N)
  mu[[1]] = cum[[1]]    # t(as.vector())

  for(n in 2:N){
    PA<-Partition_Type_All(n)
    el_js = PA$eL_r
    S_r_j<-PA$S_r_j
    an = 1:n
    L_szor_cum = array(0,c(d^n,n))
    L_szor_cum[,1] = cum[[n]]

    if ( n>=3 ) {
      for(r in 2:(n-1)){

        if(is.null(dim(el_js[[r]]))){
          el_jsr = t(as.matrix(el_js[[r]], c(1,length(el_js[[r]]))))### row vector
        }else{
          el_jsr = el_js[[r]]
        }
        Prod_M_k = array(0, c(d^n, dim(el_jsr)[1]))
        for(m in 1:dim(el_jsr)[1]){                  
          el = el_jsr[m,]
          j_el = an[el!=0]
          el_j = el[el!=0]
          Kr_1 = vector(mode = "list", length = length(j_el))
          for(k in 1:length(j_el)){
            Mc1 <-array(cum[j_el[k]], c(el_j[k],el_j[k])) 
            Kr_1[[length(j_el)-k+1]] = .KronProd(Mc1[1,])
          }
          Kr_m = .KronProd(Kr_1)
              Prod_M_k[,m] = S_r_j[[r]][m]*as.matrix(Kr_m) #L_el%*%
        }
        L_szor_cum[,r]= apply(Prod_M_k,1,sum)
      }
    }
    McN = cum[rep(1,n)]

    L_szor_cum[,n] = as.matrix(.KronProd(McN))
    xyz<-apply(L_szor_cum,1,sum)
    mu[[n]] = indx_Symmetry(xyz,d,n)
  }
  return("Mom"= mu)
}





