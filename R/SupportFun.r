
########### Support
##
## 1. KronProd
## 2. perm

## 3. KronPower
## 4. KronProdCells
## 5. KronProdVects ??????
## 6. Hermite_N
## 7. Scaling_MultiSample
## 8. Partition_Generator Section 1.4.1 Generating all Partitions
## 9. Partition_All
## 10. Commutator_Moment_eL; 2.4.3 Moment Commutators and A.2.1 Moment Commutators
##  -Use sparse matrices


#####################

.kron2 <- function(A) kronecker(A,A)

.kron3 <- function(A) kronecker(A,kronecker(A,A))

.kron4 <- function(A) kronecker(kronecker(A,A),kronecker(A,A))

.Gkd <- function(k,d) gamma((d+k)/2)/gamma(d/2)
##################################
#KronProdList
###############################
#
# Kronecker product of a list

.KronProd<-function(Mc){
  m<-length(Mc)
  Kr<-Mc[[1]]
  if (m==1) {return("Kr"=Kr)}
  for (k in 2:m) {Kr<-kronecker(Kr,Mc[[k]])}
  return("Kr"<-Kr)
}

#########################################
#Kronecker power of a vector
#################################################
# x <- c(1:4)
.KronPower <- function(x,k){
  Xk <-  rep(list(x),k)
  return(.KronProd(Xk))
}



# is.scalar
#
# Check if an element is scalar
# @examples
# a<-c(1,2)
# is.scalar(a)
.is.scalar <- function(x) is.atomic(x) && length(x) == 1L


## arrangements::permutations(1:q)
# instead of
# perm <- function(v) {
#   n <- length(v)
#   if (n == 1) v
#   else {
#     X <- NULL
#     for (i in 1:n) X <- rbind(X, cbind(v[i],perm(v[-i]) ))
#     X
#   }
# }

##############
### L22_H4 commutator matrix for Variance_of_Esti_Kurt
.L22_H4 <- function(d){
  sz <- 1
  v1 <- 1:4
  v2 <- 5:8
  perm22_H4 <- matrix(rep(0,72*8),nrow = 72,ncol = 8)
  Cperm1  <-  arrangements::combinations(4, 2)#(v1,2);
  Cperm2  <-  arrangements::combinations(v2,2);
  for (k in 1:6){
    for (j in 1:6){
      A <- c(setdiff(1:4,Cperm1[j,]), setdiff(5:8,Cperm2[k,]))
      perm22_H4[sz,]  <-  c(A[1], A[3], A[2], A[4], Cperm1[j,], Cperm2[k,])
      perm22_H4[sz+1,]  <- c(A[1], A[4], A[2], A[3], Cperm1[j,], Cperm2[k,])
      sz <- sz+2;
    }
  }
  L22_H4  <- matrix(rep(0,d^16),nrow=d^8);
  for (k in 1:72){
    L22_H4  <- L22_H4 + matr_Commutator_Kperm(perm22_H4[k,],d )
  }
  return(L22_H4)
}

#
# #' Third d-variate T-Hermite polynomial at standardized vector x
# #' @param x multivariate data of size d
.Hermite_Third<-function(x){
  d=length(x)
  H3<-Hermite_Poly_HN_Multi(x,3,diag(d))[[3]]
  return(as.vector(H3))
}

# Fourth d-variate T-Hermite polynomial at standardized vector x
# @param x multivariate data of size d
.Hermite_Fourth<-function(x){
  d=length(x)
  H4<-Hermite_Poly_HN_Multi(x,4,diag(d))[[4]]
  return(as.vector(H4))
}

#Hermite_N<-function(x,N){
#  d=length(x)
#  HN<-Hermite_Poly_HN_Multi(x,N,diag(d))[[N]]
#  return(as.vector(HN))
#}


# #' Generates all set-partition matrices for N elements
# #'
# #' @param N The (integer) number of elements to be partitioned
# #' @return U a matrix containing all the partion matrices
# #' @return S_N_r the number of partitions of size r including r=1 blocks,
# #' r=2 blocks etc.
# #' @examples
# #' Partition_Generator(3)
# #' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
# #' Springer 2021.  Section 1.4.1

.Partition_Generator<-function(N){
  U1<-matrix(c(1,1,1,0,0,1),3,2,byrow=T)
  #U1
  if (N==1) {return(list("U"=1,"S_N_r"=1))}
  if (N==2) {return(list("U"=U1, "S_N_r"=c(1,1)))}
  S_N_r <- c(1,1)
  szam2 <- NULL
  U2 <- NULL
  szam2[1] <- 1
  sep <- cumsum(S_N_r*(1:2))

  for (n in 3:N) {
    U2 <- NULL
    U2<-rbind(rep(1,n))
    U2<-rbind(U2,c(U1[1,],0))
    U2<-rbind(U2,c(rep(0,n-1),1))
    #U2
    for (m in 2:(n-1)){
      szam2[m] <-  S_N_r[m-1]+m*S_N_r[m]
      bent<-rep(0,m) #use as column vector
      for (k in 1:S_N_r[m]){
        sor <- sep[m-1]+(k-1)*m+1
        # inclusive partitions
        for (b in 1:m){
          meret<-dim(U2)
          bent[b]<-1
          U2<-rbind(U2,cbind(U1[sor:(sor+m-1),],bent))
          bent[b]<-0
        }
      }
      #Exclusive partitions
      for (k in 1:S_N_r[m]){
        sor<-sep[m-1]+(k-1)*m+1
        meret<-dim(U2)
        U2<-rbind(U2,cbind(U1[sor:(sor+m-1),],bent))
        U2<-rbind(U2,c(rep(0,n-1),1))
      }
    }
    S_N_r <- c(szam2,1)
    sep <- cumsum(S_N_r*(1:n))
    U1<-U2
    colnames(U1)<-NULL
    U2<-NULL
  }
  return(list("U"=U1,"S_N_r"=S_N_r))
}

########################
### Partition_All
# Generates all disjoint set-partition matrices for N elements
#
# #' @inheritParams Partition_Generator
# @param N The (integer) number of elements to be partitioned
# @return U the list of partition matrices  of size r.
# @return S_N_r the number of partitions of size r including r=1 blocks,
# r=2 blocks etc.
#
# @examples
# PA<-Partition_All(4)
# PA$U
# PA$S_N_r
# @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
# Springer 2021.  Section 1.4.1 (see Table 1.3, p.23, with partitions for n=4)
#
# @family Partitions
# @export
.Partition_All <- function(N){

 # if (N==8) {return(PA8)}
 # if (N==9) {return(PA9)}
 # if (N==10) {return(PA10)}

  PG<-.Partition_Generator(N)
  U1<-PG$U
  S_N_r<-PG$S_N_r
  sep<-cumsum(S_N_r*(1:N))
  U<-list()
  U[[1]]<-rep(1,N) ### first output
  if (N==1) {return(list("U"=U, "S_N_r"=1))}
  M <- 1
  for (m in 2:N) {
    U2<-U1[(sep[m-1]+1):sep[m],]
    for (k in 1:S_N_r[m]){
      M<- M+1
      sor<-(k-1)*m+1
      K<-U2[sor:(sor+m-1),]
      I<-order(apply(K,1,sum),decreasing=T)
      U[[M]]<-K[I,] # output
    }
  }
  return(list("U"=U,"S_N_r"=S_N_r))
}



# Commutator matrix  corresponding to a given  type.
#
# With an option for sparse matrices
#
# @param el_rm type of a partition
# @param d dimension
# @param useSparse T or F
# @return A commutator matrix
# @examples
# n=4;  r=2 ;  m=1 ;  d=2
# PTA<-Partition_Type_All(n)
# el_r<-PTA$eL_r[[2]][1,]
# el_r is a given type (always a vector)
# MC<- Commutator_Moment_eL(el_r,d)
# MC
# @references Gy., Terdik, Multivariate statistical methods - going beyond the linear,
# Springer 2021, Section 2.4.3, p.100, Sect. A.2.1, p. 353.,
# Corollary 2.6., p.95

##############################
### Partition_Type_eL_Location

# Find the  parameters  of a given type (eL) of partition
#
# @param eL A partition of a given type
# @return r The partition's size
# @return mm The row number in the matrix of types of
# that partition
#
#
# @examples
# eL<-c(3,1,1,0,0,0,0,0)
# Loc <-  Partition_Type_eL_Location(eL)
# # type eL is in the Loc[1]th i.e. 5th cell Loc[1]th i.e. 2nd row of types
# generated by function Partition_Type_All in $eL_r, N = 3*1+1*2+1*3=8 as follows
# PAs <-  Partition_Type_All(8)
# PAs$eL_r[[5]][2,]
# @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
# Springer 2021.  Section 1.4
#
# @family Partitions
# @export
.Partition_Type_eL_Location <- function(eL){
  N <- length(eL)
  if (sum(eL*c(1:N)) != N ){
    stop("eL is not a valid partition type")}
  eL_j<-Partition_Type_All(N)$eL_r
  loc_type_el=c(0,0)
  for (r in 1:N){
    if  (is.vector(eL_j[[r]])){
      if (prod((eL==eL_j[[r]]))){
        loc_type_el<- c(r,1)
      }
    }
    else {
      eL_jt<-eL_j[[r]]
      for (mm in 1:dim(eL_jt)[[1]]){
        if (prod((eL==eL_jt[mm,]))){
          loc_type_el<- c(r,mm)
        }
      }
    }
  }
  return(loc_type_el)
}





############################ new one!!!!!!!!!!!!!!!!
.Commutator_Moment_eL<-function(el_rm,d,useSparse=FALSE) {
  N<-length(el_rm)
  PTB<-Partition_Type_All(N)
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
    L_eL<-L_eL+matr_Commutator_Kperm(perm_Urk1,d)
    return("LeL"=L_eL)
    }
    else {
      if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
      if (.is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
      else {u_ind<-l_ind+S_m_j[[r]][m]-1}
      perm_Urk1<-matrix(0,S_m_j[[r]][m],N)
      sz<- 1
      for (k in l_ind:u_ind){
        perm_Urk1[sz,]<-Partition_2Perm(part_class[[k]])
        sz<-sz+1
      }
    }
    L_eL<-rep(0,d^N)
    for (ss in 1:dim(perm_Urk1)[1]) {
      L_eL<-L_eL + matr_Commutator_Kperm(perm_Urk1[ss,],d)
    }
  }

  if (useSparse==TRUE){
    if (r==1) {perm_Urk1<- 1:N
    L_eL<-rep(0,d^N)
    L_eL<-L_eL+matr_Commutator_Kperm(perm_Urk1,d,useSparse = TRUE)
    return("LeL"=L_eL)
    }
    else {
      if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
      if (.is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
      else {u_ind<-l_ind+S_m_j[[r]][m]-1}
      perm_Urk1<-Matrix::Matrix(0,S_m_j[[r]][m],N,sparse=TRUE)
      sz<- 1
      for (k in l_ind:u_ind){
        perm_Urk1[sz,]<-Partition_2Perm(part_class[[k]])
        sz<-sz+1
      }
    }
    L_eL<-rep(0,d^N)
    for (ss in 1:dim(perm_Urk1)[1]) {
      L_eL<-L_eL + matr_Commutator_Kperm(perm_Urk1[ss,],d,useSparse = TRUE)
    }
  }


  return("L_eL"= L_eL)
}

